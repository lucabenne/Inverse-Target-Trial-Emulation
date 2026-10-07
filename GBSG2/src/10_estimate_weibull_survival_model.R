# Fast Bayesian Weibull AFT fit using the precompiled survHEhmc Stan model.
# We bypass survHE's slow post-processing/information-criteria calculations only.
# Parameterisation: T_i ~ Weibull(shape=alpha, scale=exp(X_i beta)).

fit_weibull_hmc_fast <- function(data, formula, seed=NULL){
  for(pkg in c("survHEhmc","rstan","survHE")){
    if(!requireNamespace(pkg,quietly=TRUE)) stop("Package '",pkg,"' is required.")
  }
  ns_hmc <- asNamespace("survHEhmc")
  make_data_stan <- get("make_data_stan",envir=ns_hmc)
  stanmodels <- get("stanmodels",envir=ns_hmc)
  exArgs <- list(formula=formula,data=data)
  data.stan <- make_data_stan(formula,data,"wei",exArgs)

  if("WeibullAF" %in% names(stanmodels)){
    model_key <- "WeibullAF"
  } else {
    cand <- grep("Weibull.*AF|WeibullAF|Weibull",names(stanmodels),value=TRUE,ignore.case=TRUE)
    cand <- setdiff(cand,grep("PH",cand,value=TRUE,ignore.case=TRUE))
    if(!length(cand)) stop("No precompiled Weibull AFT model found in survHEhmc.")
    model_key <- cand[1]
  }
  dso <- stanmodels[[model_key]]
  if(is.null(seed)) seed <- sample.int(.Machine$integer.max-1L,1L)

  cat("      [direct survHEhmc Weibull] Stan sampling START\n")
  tic <- proc.time()
  fit <- rstan::sampling(
    dso,data=data.stan,chains=2,iter=2000,warmup=1000,thin=1,
    seed=seed,pars=c("beta","alpha"),include=TRUE,cores=1,
    init="random",refresh=200
  )
  cat("      [direct survHEhmc Weibull] Stan sampling DONE; elapsed =",
      round((proc.time()-tic)[3],3),"sec\n")

  ext <- rstan::extract(fit,pars=c("beta","alpha"),permuted=TRUE)
  beta_draws <- ext$beta
  if(is.null(dim(beta_draws))) beta_draws <- matrix(beta_draws,ncol=1)
  alpha_draws <- as.numeric(ext$alpha)
  if(any(!is.finite(alpha_draws)) || any(alpha_draws<=0)) stop("Invalid Weibull alpha draws.")

  list(
    beta=colMeans(beta_draws), alpha=mean(alpha_draws),
    beta_draws=beta_draws, alpha_draws=alpha_draws,
    model=fit, model_key=model_key
  )
}

# Mean latent event-time difference at a common mean covariate profile.
# The same helper works for treatment-only and covariate-adjusted formulae.
weibull_mean_survival_difference <- function(data, formula, beta, alpha){
  if(!"treat" %in% names(data)) stop("Data must contain 'treat'.")
  ref <- data[rep(1L,2L),,drop=FALSE]
  for(nm in names(ref)){
    if(nm=="treat") next
    if(is.numeric(data[[nm]]) || is.integer(data[[nm]])){
      ref[[nm]] <- rep(mean(data[[nm]],na.rm=TRUE),2L)
    }
  }
  ref$treat <- c(0,1)
  tt <- stats::delete.response(stats::terms(formula))
  X <- stats::model.matrix(tt,data=ref)
  if(ncol(X)!=length(beta)){
    stop("Reference design matrix has ",ncol(X)," columns but beta has ",length(beta)," coefficients.")
  }
  eta <- as.numeric(X %*% beta)
  mf <- gamma(1+1/alpha)
  exp(eta[2])*mf - exp(eta[1])*mf
}

# Compatibility with the original architecture.
# msd=0 returns positions [[1]]=event beta, [[2]]=censor beta, [[3]]=MSD,
# [[4]]=event alpha, [[5]]=censor alpha.
# msd=1 returns [[1]]=event beta, [[2]]=MSD, [[3]]=event alpha.
estimate_params <- function(data, formula, distrib="weibull", msd=0, seed=NULL){
  if(!tolower(distrib) %in% c("weibull","wei")) stop("Only Weibull AFT is supported in the revised GBSG2 pipeline.")
  if(is.null(seed)) seed <- sample.int(.Machine$integer.max-2L,1L)

  cat("\n      [estimate_params] Weibull event fit START\n")
  ev <- fit_weibull_hmc_fast(data,formula,seed)
  effect <- weibull_mean_survival_difference(data,formula,ev$beta,ev$alpha)
  cat("      [estimate_params] Weibull event fit DONE\n")

  if(msd==1){
    return(list(beta_event=ev$beta,mean_surv_diff=effect,alpha_event=ev$alpha))
  }

  cens <- data
  cens$event <- 1-cens$event
  cat("\n      [estimate_params] Weibull censoring fit START\n")
  ce <- fit_weibull_hmc_fast(cens,formula,seed+1L)
  cat("      [estimate_params] Weibull censoring fit DONE\n")

  list(
    beta_event=ev$beta,
    beta_cens=ce$beta,
    mean_surv_diff=effect,
    alpha_event=ev$alpha,
    alpha_cens=ce$alpha
  )
}

fit_weibull_event_and_censoring <- function(data,formula,seed=100001L){
  fit <- estimate_params(data,formula,"weibull",msd=0,seed=seed)
  list(
    event=list(beta=fit$beta_event,alpha=fit$alpha_event),
    censor=list(beta=fit$beta_cens,alpha=fit$alpha_cens),
    mean_surv_diff=fit$mean_surv_diff
  )
}

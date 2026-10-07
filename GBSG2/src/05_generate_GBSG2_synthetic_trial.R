# Baseline synthetic-RCT generator using the frozen v3 sequential covariate DGP
# and the validated Weibull AFT event/censoring DGP.
#
# IMPORTANT FOR MONTE CARLO REPLICATIONS
# --------------------------------------
# BG is fitted ONCE to the source RCT.  This file can then generate a new
# synthetic RCT using either:
#   * posterior_mode = "paired"  : legacy validated implementation; posterior
#                                  draw i is paired with synthetic subject i.
#   * posterior_mode = "dataset" : draw ONE joint posterior parameter vector
#                                  from BG and use it for the whole synthetic
#                                  dataset.  This is used by the revised
#                                  confounding Monte Carlo analysis, so each
#                                  replicate propagates DGP parameter uncertainty
#                                  without refitting JAGS.
#
# Event/censoring Weibull parameters remain the source posterior means used in
# the validated v5 pipeline.  Only the covariate-DGP posterior draw varies here.

source_if_needed <- function(file, fn){
  if(!exists(fn, mode="function")) source(file)
}

GBSG2_BG_n_draws <- function(BG){
  post <- BG$BUGSoutput$sims.list
  if(is.null(post$theta0)) stop("BG does not contain theta0 posterior draws.")
  n <- length(post$theta0)
  if(n < 1L) stop("BG contains no posterior draws.")
  n
}

resolve_GBSG2_BG_indices <- function(BG, nsim, posterior_mode=c("paired","dataset"),
                                     posterior_draw_index=NULL){
  posterior_mode <- match.arg(posterior_mode)
  npost <- GBSG2_BG_n_draws(BG)

  if(posterior_mode == "paired"){
    if(!is.null(posterior_draw_index)){
      stop("posterior_draw_index is only used with posterior_mode='dataset'.")
    }
    if(npost < nsim){
      stop("Not enough JAGS posterior draws for posterior_mode='paired': need ",
           nsim, ", have ", npost, ".")
    }
    return(seq_len(nsim))
  }

  # Dataset-level posterior predictive draw: one joint parameter vector is used
  # for every individual in the replicated synthetic RCT.
  if(is.null(posterior_draw_index)){
    posterior_draw_index <- sample.int(npost, 1L)
  }
  posterior_draw_index <- as.integer(posterior_draw_index)
  if(length(posterior_draw_index) != 1L || !is.finite(posterior_draw_index) ||
     posterior_draw_index < 1L || posterior_draw_index > npost){
    stop("Invalid posterior_draw_index: must be one integer between 1 and ", npost, ".")
  }
  rep.int(posterior_draw_index, nsim)
}


generate_GBSG2_covariates <- function(data, BG=NULL, nsim=1000L,
                                       posterior_mode=c("paired","dataset"),
                                       posterior_draw_index=NULL){
  if(!requireNamespace("truncnorm",quietly=TRUE)) stop("Package 'truncnorm' is required.")
  source_if_needed("03_fit_GBSG2_covariate_DGP.R","generate_BG")
  if(is.null(BG)) BG <- generate_BG(data)

  posterior_mode <- match.arg(posterior_mode)
  post <- BG$BUGSoutput$sims.list
  ii <- resolve_GBSG2_BG_indices(
    BG=BG,
    nsim=nsim,
    posterior_mode=posterior_mode,
    posterior_draw_index=posterior_draw_index
  )

  source_means <- c(
    age=mean(data$age), menostat=mean(data$menostat), tsize=mean(data$tsize),
    tgrade=mean(data$tgrade), pnodes=mean(data$pnodes),
    sqrt_progrec=mean(sqrt(data$progrec))
  )

  age <- truncnorm::rtruncnorm(nsim,20,Inf,post$theta0[ii],post$sigma0[ii])

  menostat <- rbinom(nsim,1,pnorm(
    post$theta10[ii] + post$theta11[ii]*(age-source_means[["age"]])
  ))

  tsize <- truncnorm::rtruncnorm(
    nsim,0,Inf,
    post$theta2[ii] +
      post$theta21[ii]*(age-source_means[["age"]]) +
      post$theta22[ii]*(menostat-source_means[["menostat"]]),
    post$sigma2[ii]
  )

  eta_grade <-
    post$theta31[ii]*(age-source_means[["age"]]) +
    post$theta32[ii]*(menostat-source_means[["menostat"]]) +
    post$theta33[ii]*(tsize-source_means[["tsize"]])
  cdf1 <- plogis(post$grade_cut1[ii]-eta_grade)
  cdf2 <- plogis(post$grade_cut2[ii]-eta_grade)
  p1 <- pmax(cdf1,0)
  p2 <- pmax(cdf2-cdf1,0)
  p3 <- pmax(1-cdf2,0)
  ps <- p1+p2+p3
  p1 <- p1/ps; p2 <- p2/ps; p3 <- p3/ps
  tgrade <- vapply(seq_len(nsim),function(k){
    sample.int(3L,1L,prob=c(p1[k],p2[k],p3[k]))
  },integer(1))

  pnodes <- truncnorm::rtruncnorm(
    nsim,0,Inf,
    post$theta40[ii] +
      post$theta41[ii]*(tsize-source_means[["tsize"]]) +
      post$theta42[ii]*(age-source_means[["age"]]) +
      post$theta43[ii]*(menostat-source_means[["menostat"]]) +
      post$theta44[ii]*(tgrade-source_means[["tgrade"]]),
    post$sigma4[ii]
  )

  sqrt_progrec <- truncnorm::rtruncnorm(
    nsim,0,Inf,
    post$theta5[ii] +
      post$theta51[ii]*(age-source_means[["age"]]) +
      post$theta52[ii]*(menostat-source_means[["menostat"]]) +
      post$theta53[ii]*(tsize-source_means[["tsize"]]) +
      post$theta54[ii]*(tgrade-source_means[["tgrade"]]) +
      post$theta55[ii]*(pnodes-source_means[["pnodes"]]),
    post$sigma5[ii]
  )
  progrec <- sqrt_progrec^2

  sqrt_estrec <- truncnorm::rtruncnorm(
    nsim,0,Inf,
    post$theta6[ii] +
      post$theta61[ii]*(age-source_means[["age"]]) +
      post$theta62[ii]*(menostat-source_means[["menostat"]]) +
      post$theta63[ii]*(tsize-source_means[["tsize"]]) +
      post$theta64[ii]*(tgrade-source_means[["tgrade"]]) +
      post$theta65[ii]*(pnodes-source_means[["pnodes"]]) +
      post$theta66[ii]*(sqrt_progrec-source_means[["sqrt_progrec"]]),
    post$sigma6[ii]
  )
  estrec <- sqrt_estrec^2

  # Preserve the integer-valued GBSG2 implementation used in the original code.
  out <- round(data.frame(age,menostat,tsize,tgrade,pnodes,progrec,estrec))
  if(any(!out$tgrade %in% 1:3)) stop("Generated tgrade outside {1,2,3}.")

  attr(out,"posterior_mode") <- posterior_mode
  if(posterior_mode == "dataset") attr(out,"BG_draw_index") <- ii[1]
  out
}


generate_trial <- function(data, beta_weib, beta_weib_cens,
                           alpha_event, alpha_cens, scaling,
                           BG=NULL, nsim=1000L,
                           posterior_mode=c("paired","dataset"),
                           posterior_draw_index=NULL){
  if(!exists("GBSG2_draw_weibull_time",mode="function")) source("00_GBSG2_helpers.R")

  posterior_mode <- match.arg(posterior_mode)
  covariates <- generate_GBSG2_covariates(
    data=data,
    BG=BG,
    nsim=nsim,
    posterior_mode=posterior_mode,
    posterior_draw_index=posterior_draw_index
  )
  used_bg_draw <- attr(covariates,"BG_draw_index")

  treat <- rbinom(nsim,1,0.5)
  sim_data <- data.frame(treat=treat,covariates)

  event_time <- GBSG2_draw_weibull_time(sim_data,beta_weib,alpha_event,scaling)
  censor_time <- GBSG2_draw_weibull_time(sim_data,beta_weib_cens,alpha_cens,scaling)

  sim_data$time <- pmin(event_time,censor_time)
  sim_data$event <- as.numeric(event_time<censor_time)

  attr(sim_data,"posterior_mode") <- posterior_mode
  if(posterior_mode == "dataset") attr(sim_data,"BG_draw_index") <- used_bg_draw
  sim_data
}

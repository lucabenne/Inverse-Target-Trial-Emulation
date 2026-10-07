# Additional Sampling bias mechanism updated to the frozen v3 covariate DGP
# and Weibull AFT outcome/censoring model.
#
# The bias mechanism itself is unchanged: select a distant treated/control seed
# pair and append locally generated subjects around those profiles.
#
# For Monte Carlo replication BG is still fitted only once.  The caller may pass
# BG_draw_index to use one joint posterior DGP draw for the whole replicate.
# If BG_draw_index is NULL, posterior means are used (legacy v5 behaviour).

newdatamahal_distr_surv2 <- function(data, bern, n, q,
                                     beta_weib, beta_weib_cens,
                                     alpha_event, alpha_cens, scaling,
                                     Ma, BG=NULL, data_prec=NULL,
                                     BG_draw_index=NULL){
  if(!requireNamespace("truncnorm",quietly=TRUE)) stop("Package 'truncnorm' is required.")
  if(!exists("GBSG2_draw_weibull_time",mode="function")) source("00_GBSG2_helpers.R")
  if(!exists("generate_BG",mode="function")) source("03_fit_GBSG2_covariate_DGP.R")
  if(is.null(BG)) BG <- generate_BG(data)
  if(is.null(data_prec)) data_prec <- data

  n <- as.integer(n)
  if(length(n)!=1L || n<1L) stop("n must be one positive integer.")
  if(length(bern)<n) stop("bern must contain at least n treatment draws.")

  cov_names <- c("age","menostat","tsize","tgrade","pnodes","progrec","estrec")
  ind0 <- which(data$treat==0)
  ind1 <- which(data$treat==1)
  post <- BG$BUGSoutput$sims.list

  if(is.null(BG_draw_index)){
    # Legacy behaviour: local conditional slopes/cutpoints are posterior means.
    dgp <- vapply(post,mean,numeric(1))
  } else {
    BG_draw_index <- as.integer(BG_draw_index)
    npost <- length(post$theta0)
    if(length(BG_draw_index)!=1L || !is.finite(BG_draw_index) ||
       BG_draw_index<1L || BG_draw_index>npost){
      stop("BG_draw_index must be one integer between 1 and ",npost,".")
    }
    # One JOINT posterior vector for the whole Monte Carlo replicate.
    dgp <- vapply(post,function(x) as.numeric(x[BG_draw_index]),numeric(1))
  }

  source_means <- c(
    age=mean(data$age), menostat=mean(data$menostat), tsize=mean(data$tsize),
    tgrade=mean(data$tgrade), pnodes=mean(data$pnodes),
    sqrt_progrec=mean(sqrt(data$progrec)), sqrt_estrec=mean(sqrt(data$estrec))
  )

  local_sd <- function(seed,center,minimum){
    max(abs(seed-center)/4,minimum)
  }

  grade_draw <- function(age,menostat,tsize){
    eta <- dgp[["theta31"]]*(age-source_means[["age"]]) +
      dgp[["theta32"]]*(menostat-source_means[["menostat"]]) +
      dgp[["theta33"]]*(tsize-source_means[["tsize"]])
    c1 <- plogis(dgp[["grade_cut1"]]-eta)
    c2 <- plogis(dgp[["grade_cut2"]]-eta)
    pr <- pmax(c(c1,c2-c1,1-c2),0)
    pr <- pr/sum(pr)
    sample.int(3L,1L,prob=pr)
  }

  generate_near_seed <- function(seed,treat_value){
    out <- as.list(seed)

    out$age <- truncnorm::rtruncnorm(
      1,20,Inf,seed$age,local_sd(seed$age,source_means[["age"]],0.5)
    )

    # Preserve the original local anchoring of menopausal status while using
    # the fitted age slope from the selected full-sequential posterior draw.
    latent_seed <- if(seed$menostat==1) qnorm(0.75) else qnorm(0.25)
    p_meno <- pnorm(latent_seed + dgp[["theta11"]]*(out$age-seed$age))
    out$menostat <- rbinom(1,1,p_meno)

    mu_tsize <- seed$tsize +
      dgp[["theta21"]]*(out$age-seed$age) +
      dgp[["theta22"]]*(out$menostat-seed$menostat)
    out$tsize <- truncnorm::rtruncnorm(
      1,0,Inf,mu_tsize,local_sd(seed$tsize,source_means[["tsize"]],0.5)
    )

    out$tgrade <- grade_draw(out$age,out$menostat,out$tsize)

    mu_pnodes <- seed$pnodes +
      dgp[["theta41"]]*(out$tsize-seed$tsize) +
      dgp[["theta42"]]*(out$age-seed$age) +
      dgp[["theta43"]]*(out$menostat-seed$menostat) +
      dgp[["theta44"]]*(out$tgrade-seed$tgrade)
    out$pnodes <- truncnorm::rtruncnorm(
      1,0,Inf,mu_pnodes,local_sd(seed$pnodes,source_means[["pnodes"]],0.25)
    )

    seed_sp <- sqrt(max(seed$progrec,0))
    mu_sp <- seed_sp +
      dgp[["theta51"]]*(out$age-seed$age) +
      dgp[["theta52"]]*(out$menostat-seed$menostat) +
      dgp[["theta53"]]*(out$tsize-seed$tsize) +
      dgp[["theta54"]]*(out$tgrade-seed$tgrade) +
      dgp[["theta55"]]*(out$pnodes-seed$pnodes)
    sp <- truncnorm::rtruncnorm(
      1,0,Inf,mu_sp,local_sd(seed_sp,source_means[["sqrt_progrec"]],0.1)
    )
    out$progrec <- sp^2

    seed_se <- sqrt(max(seed$estrec,0))
    mu_se <- seed_se +
      dgp[["theta61"]]*(out$age-seed$age) +
      dgp[["theta62"]]*(out$menostat-seed$menostat) +
      dgp[["theta63"]]*(out$tsize-seed$tsize) +
      dgp[["theta64"]]*(out$tgrade-seed$tgrade) +
      dgp[["theta65"]]*(out$pnodes-seed$pnodes) +
      dgp[["theta66"]]*(sp-seed_sp)
    se <- truncnorm::rtruncnorm(
      1,0,Inf,mu_se,local_sd(seed_se,source_means[["sqrt_estrec"]],0.1)
    )
    out$estrec <- se^2

    row_cov <- data.frame(
      treat=treat_value,
      age=out$age,menostat=out$menostat,tsize=out$tsize,tgrade=out$tgrade,
      pnodes=out$pnodes,progrec=out$progrec,estrec=out$estrec
    )
    te <- GBSG2_draw_weibull_time(row_cov,beta_weib,alpha_event,scaling)
    tc <- GBSG2_draw_weibull_time(row_cov,beta_weib_cens,alpha_cens,scaling)
    row_cov$time <- min(te,tc)
    row_cov$event <- as.numeric(te<tc)
    row_cov
  }

  for(j in seq_along(q)){
    Ma_ordered <- sort(Ma)
    qq <- min(max(q[j],0),1)
    pos <- max(1L,min(length(Ma_ordered),round(length(Ma_ordered)*qq)))
    threshold <- Ma_ordered[pos]
    selected_pair <- which(Ma==threshold,arr.ind=TRUE)[1,]
    seed_control <- data[ind0[selected_pair[1]],,drop=FALSE]
    seed_treated <- data[ind1[selected_pair[2]],,drop=FALSE]

    for(i in seq_len(n)){
      if(bern[i]==1){
        nr <- generate_near_seed(seed_treated,1)
      } else {
        nr <- generate_near_seed(seed_control,0)
      }
      nr <- nr[,names(data),drop=FALSE]
      data_prec <- rbind(data_prec,nr)
    }
  }

  # Preserve the integer-valued implementation used in the original scripts.
  data_prec <- round(data_prec)
  if(any(!data_prec$tgrade %in% 1:3)) stop("Additional Sampling generated invalid tgrade.")
  data_prec
}

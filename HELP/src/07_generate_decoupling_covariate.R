decoupling_covariate = function(data,N,p,x,BG){

  library(readr)
  library(MatchIt)

  if (!exists("generate_HELP_synthetic_trial")) {
    source("03a_generate_HELP_synthetic_trial.R")
  }

  # Prognostic score coefficients are estimated once from the source HELP RCT.
  std_vars = c("dayslink","daysdrink","daysanysub","pcs","mcs","cesd","sexrisk","drugrisk")
  source_means = sapply(data[std_vars],mean,na.rm=TRUE)
  source_sds = sapply(data[std_vars],sd,na.rm=TRUE)

  data_std = data
  for (v in std_vars) {
    data_std[[v]] = (data[[v]]-source_means[v])/source_sds[v]
  }

  m_std = lm(
    dayslink ~ treat + female + daysdrink + daysanysub + pcs + mcs + cesd + sexrisk + drugrisk,
    data = data_std
  )

  beta_female = unname(coef(m_std)["female"])
  beta_x = unname(coef(m_std)[x])

  make_one = function(N_target) {

    # Fresh randomized synthetic RCT for this replicate. The DGP itself was fit
    # once; baseline-covariate parameters use one joint posterior draw here,
    # while outcome-model parameters stay fixed at posterior means.
    data1 = generate_HELP_synthetic_trial(
      data = data,
      BG = BG,
      nsim = 1000L,
      treat_prob = p,
      posterior_mode = "dataset",
      outcome_posterior_mode = "mean"
    )

    target_rct = attr(data1, "rct_target")

    score = beta_female*data1$female +
      beta_x*((data1[[x]]-source_means[x])/source_sds[x])

    expected_removals = function(alpha,gamma){
      eta = alpha + gamma*score
      keep_treat = exp(pmin(0,eta[data1$treat==1]))
      keep_control = exp(pmin(0,-eta[data1$treat==0]))
      c(
        treated = sum(1-keep_treat),
        control = sum(1-keep_control)
      )
    }

    find_alpha = function(gamma){
      f = function(alpha){
        er = expected_removals(alpha,gamma)
        er["treated"]-er["control"]
      }

      lo = -5
      hi = 5
      while (sign(f(lo)) == sign(f(hi))) {
        lo = lo*2
        hi = hi*2
      }
      uniroot(f,interval=c(lo,hi),tol=1e-10)$root
    }

    calibrate = function(target){
      gfun = function(gamma){
        alpha = find_alpha(gamma)
        mean(expected_removals(alpha,gamma))-target
      }

      gamma_lo = 1e-8
      gamma_hi = 1
      while (gfun(gamma_hi) < 0) gamma_hi = gamma_hi*2

      gamma = uniroot(gfun,interval=c(gamma_lo,gamma_hi),tol=1e-10)$root
      alpha = find_alpha(gamma)
      c(alpha=alpha,gamma=gamma)
    }

    cal = calibrate(N_target)
    alpha = cal["alpha"]
    gamma = cal["gamma"]

    eta = alpha + gamma*score
    p_keep = ifelse(
      data1$treat==1,
      exp(pmin(0,eta)),
      exp(pmin(0,-eta))
    )

    keep = rbinom(nrow(data1),1,p_keep)==1
    out = data1[keep,,drop=FALSE]
    rownames(out) = NULL
    attr(out, "rct_target") = target_rct
    out
  }

  if (length(N) == 1L) return(make_one(N))

  list_data = vector("list",length(N))
  for (k in seq_along(N)) {
    if (k == 1L || k %% 25L == 0L || k == length(N)) {
      cat("Decoupling Covariate replicate",k,"of",length(N),"\n")
    }
    list_data[[k]] = make_one(N[k])
  }
  list_data
}

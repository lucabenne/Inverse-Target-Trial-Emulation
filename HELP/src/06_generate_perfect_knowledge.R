library(readr)
library(truncnorm)

source("../src/10_estimate_IPW.R")
source("../src/08_measure_mahalanobis_imbalance.R")
source("../src/12_estimate_AIPW.R")
if (!exists("HELP_generate_covariates")) source("../src/03a_generate_HELP_synthetic_trial.R")

observational_help2 = function(data, n, BG, gamma = 1.25){

  list_data = list()

  # Prognostic score estimated once from the source HELP RCT.
  # female and daysanysub define the known treatment-allocation mechanism.
  std_vars = c("dayslink","daysdrink","daysanysub","pcs","mcs","cesd","sexrisk","drugrisk")
  source_means = sapply(data[std_vars], mean, na.rm = TRUE)
  source_sds = sapply(data[std_vars], sd, na.rm = TRUE)

  data_std = data
  for (v in std_vars) {
    data_std[[v]] = (data[[v]] - source_means[v])/source_sds[v]
  }

  m_std = lm(
    dayslink ~ treat + female + daysdrink + daysanysub + pcs + mcs + cesd + sexrisk + drugrisk,
    data = data_std
  )

  beta_female = unname(coef(m_std)["female"])
  beta_daysanysub = unname(coef(m_std)["daysanysub"])

  source_score_raw = beta_female*data$female +
    beta_daysanysub*((data$daysanysub-source_means["daysanysub"])/source_sds["daysanysub"])

  source_score_mean = mean(source_score_raw)
  source_score_sd = sd(source_score_raw)
  source_score = (source_score_raw-source_score_mean)/source_score_sd

  # Confounding strength is supplied explicitly by the analysis script.
  # The moderate scenario uses gamma = 1.25.

  # Calibrate the intercept so expected treatment prevalence is approximately .50
  # in the source covariate distribution.
  f_alpha = function(alpha){
    mean(plogis(alpha + gamma*source_score)) - 0.50
  }
  alpha = uniroot(f_alpha, interval = c(-20,20), tol = 1e-12)$root

  for (j in 1:length(n)) {

    # One joint posterior draw for the baseline-covariate DGP in this replicate.
    # Outcome-model parameters are fixed at posterior means, as in GBSG2.
    posterior_draw = HELP_get_posterior_draw(BG, "dataset")

    covariates = HELP_generate_covariates(
      data = data,
      BG = BG,
      nsim = n[j],
      posterior_mode = "dataset",
      posterior_draw = posterior_draw
    )

    # Perfect Knowledge: treatment comes directly from the known logistic PS.
    score_raw = beta_female*covariates$female +
      beta_daysanysub*((covariates$daysanysub-source_means["daysanysub"])/source_sds["daysanysub"])
    score = (score_raw-source_score_mean)/source_score_sd

    p = plogis(alpha + gamma*score)
    treat = rbinom(n[j],1,p)

    dayslink = HELP_generate_outcome(
      covariates = covariates,
      treat = treat,
      BG = BG,
      posterior_mode = "mean",
      posterior_draw = NULL,
      round_outcome = FALSE
    )

    data_sim = data.frame(
      treat = treat,
      female = covariates$female,
      dayslink = dayslink,
      daysdrink = covariates$daysdrink,
      daysanysub = covariates$daysanysub,
      pcs = covariates$pcs,
      mcs = covariates$mcs,
      cesd = covariates$cesd,
      sexrisk = covariates$sexrisk,
      drugrisk = covariates$drugrisk
    )

    # Replicate-specific target: randomized synthetic RCT on the SAME generated
    # covariate population, before making treatment depend on X.  This auxiliary
    # benchmark is generated without advancing the main Monte Carlo RNG stream.
    rng_state = .Random.seed
    treat_rct = rbinom(n[j], 1, 0.5)
    dayslink_rct = HELP_generate_outcome(
      covariates = covariates,
      treat = treat_rct,
      BG = BG,
      posterior_mode = "mean",
      posterior_draw = NULL,
      round_outcome = FALSE
    )
    rct_target = mean(dayslink_rct[treat_rct == 0]) -
      mean(dayslink_rct[treat_rct == 1])
    assign(".Random.seed", rng_state, envir = .GlobalEnv)

    attr(data_sim, "rct_target") = rct_target
    list_data = append(list_data,list(data_sim))
  }

  return(list_data)
}

library(truncnorm)

# -----------------------------------------------------------------------------
# Shared HELP synthetic-data generator.
#
# The baseline covariates are generated from the sequential conditional DGP
# fitted in 02_fit_HELP_bayesian_DGP.R / 03_HELP_DGP_model.txt.
#
# posterior_mode controls the baseline-covariate DGP. In the paper pipeline,
# each replicate uses one joint posterior draw for all covariate-model parameters.
# outcome_posterior_mode controls the outcome DGP separately. The paper pipeline
# fixes outcome parameters at posterior means so the treatment-effect target is
# constant across Monte Carlo replicates, matching the GBSG2 architecture.
# -----------------------------------------------------------------------------

HELP_get_posterior_draw <- function(BG,
                                    posterior_mode = c("dataset", "mean"),
                                    posterior_draw = NULL) {
  posterior_mode <- match.arg(posterior_mode)
  sims <- BG$BUGSoutput$sims.list
  n_draws <- length(sims$theta90)

  if (posterior_mode == "mean") return(NA_integer_)

  if (is.null(posterior_draw)) {
    posterior_draw <- sample.int(n_draws, 1L)
  }

  posterior_draw <- as.integer(posterior_draw)
  if (posterior_draw < 1L || posterior_draw > n_draws) {
    stop("posterior_draw is outside the available posterior draws.")
  }

  posterior_draw
}

HELP_get_scalar <- function(BG, name, posterior_mode, posterior_draw) {
  x <- BG$BUGSoutput$sims.list[[name]]
  if (is.null(x)) stop("Posterior parameter not found: ", name)
  if (posterior_mode == "mean") mean(x) else x[posterior_draw]
}

HELP_get_vector <- function(BG, name, posterior_mode, posterior_draw) {
  x <- BG$BUGSoutput$sims.list[[name]]
  if (is.null(x)) stop("Posterior parameter not found: ", name)

  if (is.null(dim(x))) {
    return(if (posterior_mode == "mean") mean(x) else x[posterior_draw])
  }

  if (posterior_mode == "mean") {
    as.numeric(apply(x, 2L, mean))
  } else {
    as.numeric(x[posterior_draw, ])
  }
}

HELP_standardize <- function(x, variable, BG) {
  m <- BG$HELP_scaling$means[variable]
  s <- BG$HELP_scaling$sds[variable]
  (x - m) / s
}

HELP_crude_effect <- function(d) {
  if (!all(c("treat", "dayslink") %in% names(d))) {
    stop("Dataset must contain treat and dayslink.")
  }
  if (!all(c(0, 1) %in% unique(d$treat))) {
    return(NA_real_)
  }
  mean(d$dayslink[d$treat == 0]) - mean(d$dayslink[d$treat == 1])
}

HELP_generate_covariates <- function(data,
                                     BG,
                                     nsim = 1000L,
                                     posterior_mode = c("dataset", "mean"),
                                     posterior_draw = NULL,
                                     seed = NULL) {

  posterior_mode <- match.arg(posterior_mode)
  if (!is.null(seed)) set.seed(seed)

  if (is.null(BG$HELP_scaling)) {
    stop("BG does not contain HELP_scaling. Refit the HELP DGP first.")
  }

  nsim <- as.integer(nsim)
  if (length(nsim) != 1L || nsim < 1L) stop("nsim must be a positive integer.")

  posterior_draw <- HELP_get_posterior_draw(
    BG = BG,
    posterior_mode = posterior_mode,
    posterior_draw = posterior_draw
  )

  theta0 <- HELP_get_scalar(BG, "theta0", posterior_mode, posterior_draw)
  female <- rbinom(nsim, 1, theta0)
  z_female <- HELP_standardize(female, "female", BG)

  b <- HELP_get_vector(BG, "beta_drink", posterior_mode, posterior_draw)
  mu <- b[1] + b[2] * z_female
  daysdrink <- round(rtruncnorm(
    nsim, a = 0, b = Inf, mean = mu,
    sd = HELP_get_scalar(BG, "sigma_daysdrink", posterior_mode, posterior_draw)
  ))
  z_daysdrink <- HELP_standardize(daysdrink, "daysdrink", BG)

  b <- HELP_get_vector(BG, "beta_sub", posterior_mode, posterior_draw)
  mu <- b[1] + b[2] * z_female + b[3] * z_daysdrink
  daysanysub <- round(rtruncnorm(
    nsim, a = 0, b = Inf, mean = mu,
    sd = HELP_get_scalar(BG, "sigma_daysanysub", posterior_mode, posterior_draw)
  ))
  z_daysanysub <- HELP_standardize(daysanysub, "daysanysub", BG)

  b <- HELP_get_vector(BG, "beta_pcs", posterior_mode, posterior_draw)
  mu <- b[1] + b[2] * z_female + b[3] * z_daysdrink + b[4] * z_daysanysub
  pcs <- rtruncnorm(
    nsim, a = 14, b = 75, mean = mu,
    sd = HELP_get_scalar(BG, "sigma_pcs", posterior_mode, posterior_draw)
  )
  z_pcs <- HELP_standardize(pcs, "pcs", BG)

  b <- HELP_get_vector(BG, "beta_mcs", posterior_mode, posterior_draw)
  mu <- b[1] + b[2] * z_female + b[3] * z_daysdrink +
    b[4] * z_daysanysub + b[5] * z_pcs
  mcs <- rtruncnorm(
    nsim, a = 6.6, b = 62, mean = mu,
    sd = HELP_get_scalar(BG, "sigma_mcs", posterior_mode, posterior_draw)
  )
  z_mcs <- HELP_standardize(mcs, "mcs", BG)

  b <- HELP_get_vector(BG, "beta_cesd", posterior_mode, posterior_draw)
  mu <- b[1] + b[2] * z_female + b[3] * z_daysdrink +
    b[4] * z_daysanysub + b[5] * z_pcs + b[6] * z_mcs
  cesd <- round(rtruncnorm(
    nsim, a = 0, b = 60, mean = mu,
    sd = HELP_get_scalar(BG, "sigma_cesd", posterior_mode, posterior_draw)
  ))
  cesd <- pmin(pmax(cesd, 0), 60)
  z_cesd <- HELP_standardize(cesd, "cesd", BG)

  b <- HELP_get_vector(BG, "beta_sexrisk", posterior_mode, posterior_draw)
  mu <- b[1] + b[2] * z_female + b[3] * z_daysdrink +
    b[4] * z_daysanysub + b[5] * z_pcs + b[6] * z_mcs + b[7] * z_cesd
  sexrisk <- round(rtruncnorm(
    nsim, a = 0, b = 21, mean = mu,
    sd = HELP_get_scalar(BG, "sigma_sexrisk", posterior_mode, posterior_draw)
  ))
  sexrisk <- pmin(pmax(sexrisk, 0), 21)
  z_sexrisk <- HELP_standardize(sexrisk, "sexrisk", BG)

  b <- HELP_get_vector(BG, "beta_drugrisk", posterior_mode, posterior_draw)
  mu <- b[1] + b[2] * z_female + b[3] * z_daysdrink +
    b[4] * z_daysanysub + b[5] * z_pcs + b[6] * z_mcs +
    b[7] * z_cesd + b[8] * z_sexrisk
  drugrisk <- round(rtruncnorm(
    nsim, a = 0, b = 21, mean = mu,
    sd = HELP_get_scalar(BG, "sigma_drugrisk", posterior_mode, posterior_draw)
  ))
  drugrisk <- pmin(pmax(drugrisk, 0), 21)

  out <- data.frame(
    female = female,
    daysdrink = daysdrink,
    daysanysub = daysanysub,
    pcs = pcs,
    mcs = mcs,
    cesd = cesd,
    sexrisk = sexrisk,
    drugrisk = drugrisk
  )

  attr(out, "posterior_mode") <- posterior_mode
  attr(out, "posterior_draw") <- posterior_draw
  out
}

HELP_generate_outcome <- function(covariates,
                                  treat,
                                  BG,
                                  posterior_mode = c("dataset", "mean"),
                                  posterior_draw = NULL,
                                  round_outcome = FALSE) {

  posterior_mode <- match.arg(posterior_mode)
  covariates <- as.data.frame(covariates)
  treat <- as.numeric(treat)

  if (nrow(covariates) != length(treat)) {
    stop("covariates and treat must have the same number of rows.")
  }

  if (posterior_mode == "dataset" && is.null(posterior_draw)) {
    posterior_draw <- HELP_get_posterior_draw(BG, posterior_mode)
  }

  z_female <- HELP_standardize(covariates$female, "female", BG)
  z_daysdrink <- HELP_standardize(covariates$daysdrink, "daysdrink", BG)
  z_daysanysub <- HELP_standardize(covariates$daysanysub, "daysanysub", BG)
  z_pcs <- HELP_standardize(covariates$pcs, "pcs", BG)
  z_mcs <- HELP_standardize(covariates$mcs, "mcs", BG)
  z_cesd <- HELP_standardize(covariates$cesd, "cesd", BG)
  z_sexrisk <- HELP_standardize(covariates$sexrisk, "sexrisk", BG)
  z_drugrisk <- HELP_standardize(covariates$drugrisk, "drugrisk", BG)

  mu <-
    HELP_get_scalar(BG, "theta90", posterior_mode, posterior_draw) +
    HELP_get_scalar(BG, "theta91", posterior_mode, posterior_draw) * treat +
    HELP_get_scalar(BG, "theta92", posterior_mode, posterior_draw) * z_daysdrink +
    HELP_get_scalar(BG, "theta93", posterior_mode, posterior_draw) * z_daysanysub +
    HELP_get_scalar(BG, "theta94", posterior_mode, posterior_draw) * z_female +
    HELP_get_scalar(BG, "theta95", posterior_mode, posterior_draw) * z_pcs +
    HELP_get_scalar(BG, "theta96", posterior_mode, posterior_draw) * z_mcs +
    HELP_get_scalar(BG, "theta97", posterior_mode, posterior_draw) * z_cesd +
    HELP_get_scalar(BG, "theta98", posterior_mode, posterior_draw) * z_sexrisk +
    HELP_get_scalar(BG, "theta99", posterior_mode, posterior_draw) * z_drugrisk

  y <- rtruncnorm(
    nrow(covariates), a = 0, b = Inf, mean = mu,
    sd = HELP_get_scalar(BG, "sigma_dayslink", posterior_mode, posterior_draw)
  )

  # T(0, ) implies a strictly positive outcome.  When integer days are
  # requested (Additional Sampling), ordinary round() can turn values in
  # (0, 0.5) into exactly 0, which is outside the support of the Gamma
  # models used in the deliberate outcome-misspecification analysis.
  if (round_outcome) {
    y <- pmax(1, round(y))
  } else {
    y <- pmax(y, .Machine$double.eps)
  }
  y
}

generate_HELP_synthetic_trial <- function(data,
                                          BG,
                                          nsim = 1000L,
                                          treat_prob = mean(data$treat),
                                          posterior_mode = c("dataset", "mean"),
                                          posterior_draw = NULL,
                                          outcome_posterior_mode = c("mean", "dataset"),
                                          seed = NULL) {

  posterior_mode <- match.arg(posterior_mode)
  outcome_posterior_mode <- match.arg(outcome_posterior_mode)
  if (!is.null(seed)) set.seed(seed)

  if (!is.numeric(treat_prob) || length(treat_prob) != 1L ||
      treat_prob <= 0 || treat_prob >= 1) {
    stop("treat_prob must be strictly between 0 and 1.")
  }

  posterior_draw <- HELP_get_posterior_draw(
    BG = BG,
    posterior_mode = posterior_mode,
    posterior_draw = posterior_draw
  )

  covariates <- HELP_generate_covariates(
    data = data,
    BG = BG,
    nsim = nsim,
    posterior_mode = posterior_mode,
    posterior_draw = posterior_draw
  )

  treat <- rbinom(nsim, 1, treat_prob)

  outcome_draw <- if (outcome_posterior_mode == "dataset") posterior_draw else NULL

  dayslink <- HELP_generate_outcome(
    covariates = covariates,
    treat = treat,
    BG = BG,
    posterior_mode = outcome_posterior_mode,
    posterior_draw = outcome_draw,
    round_outcome = FALSE
  )

  out <- data.frame(
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

  attr(out, "covariate_posterior_mode") <- posterior_mode
  attr(out, "covariate_posterior_draw") <- posterior_draw
  attr(out, "outcome_posterior_mode") <- outcome_posterior_mode
  attr(out, "treat_prob") <- treat_prob
  # Replicate-specific randomized-trial benchmark before any ITTE mechanism.
  attr(out, "rct_target") <- HELP_crude_effect(out)
  out
}

library(R2jags)

HELP_required_variables <- c(
  "treat", "female", "dayslink", "daysdrink", "daysanysub",
  "pcs", "mcs", "cesd", "sexrisk", "drugrisk"
)

HELP_covariates <- c(
  "female", "daysdrink", "daysanysub", "pcs", "mcs",
  "cesd", "sexrisk", "drugrisk"
)

HELP_safe_sd <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) stop("A source variable has zero/invalid SD.")
  s
}

HELP_source_scaling <- function(data) {
  vars <- c(HELP_covariates, "dayslink")
  list(
    means = sapply(data[vars], mean, na.rm = TRUE),
    sds   = sapply(data[vars], HELP_safe_sd)
  )
}

jags_data <- function(data,
                      model_file = "03_HELP_DGP_model.txt",
                      seed = 101010L,
                      n.chains = 2L,
                      n.iter = 3000L,
                      n.burnin = 1000L,
                      n.thin = 2L) {

  missing_vars <- setdiff(HELP_required_variables, names(data))
  if (length(missing_vars) > 0L) {
    stop("Missing required HELP variables: ", paste(missing_vars, collapse = ", "))
  }

  data <- as.data.frame(data[, HELP_required_variables])
  if (anyNA(data)) stop("The HELP DGP fit requires complete data for the selected variables.")

  N <- nrow(data)
  if (N < 2L) stop("Not enough observations to fit the HELP DGP.")

  scaling <- HELP_source_scaling(data)
  m <- scaling$means
  s <- scaling$sds

  # Predictors are standardized.  For each conditional response, slopes receive
  # weak Normal priors with SD equal to twice the source SD of that response.
  prior_prec_daysdrink  <- 1 / (2 * s["daysdrink"])^2
  prior_prec_daysanysub <- 1 / (2 * s["daysanysub"])^2
  prior_prec_pcs        <- 1 / (2 * s["pcs"])^2
  prior_prec_mcs        <- 1 / (2 * s["mcs"])^2
  prior_prec_cesd       <- 1 / (2 * s["cesd"])^2
  prior_prec_sexrisk    <- 1 / (2 * s["sexrisk"])^2
  prior_prec_drugrisk   <- 1 / (2 * s["drugrisk"])^2
  prior_prec_dayslink   <- 1 / (2 * s["dayslink"])^2

  datalist <- list(
    N = N,
    treat = data$treat,
    female = data$female,
    dayslink = data$dayslink,
    daysdrink = data$daysdrink,
    daysanysub = data$daysanysub,
    pcs = data$pcs,
    mcs = data$mcs,
    cesd = data$cesd,
    sexrisk = data$sexrisk,
    drugrisk = data$drugrisk,

    center_female = m["female"],
    center_daysdrink = m["daysdrink"],
    center_daysanysub = m["daysanysub"],
    center_pcs = m["pcs"],
    center_mcs = m["mcs"],
    center_cesd = m["cesd"],
    center_sexrisk = m["sexrisk"],
    center_drugrisk = m["drugrisk"],

    scale_female = s["female"],
    scale_daysdrink = s["daysdrink"],
    scale_daysanysub = s["daysanysub"],
    scale_pcs = s["pcs"],
    scale_mcs = s["mcs"],
    scale_cesd = s["cesd"],
    scale_sexrisk = s["sexrisk"],
    scale_drugrisk = s["drugrisk"],

    prior_mean_daysdrink = m["daysdrink"],
    prior_mean_daysanysub = m["daysanysub"],
    prior_mean_pcs = m["pcs"],
    prior_mean_mcs = m["mcs"],
    prior_mean_cesd = m["cesd"],
    prior_mean_sexrisk = m["sexrisk"],
    prior_mean_drugrisk = m["drugrisk"],
    prior_mean_dayslink = m["dayslink"],

    prior_prec_daysdrink = prior_prec_daysdrink,
    prior_prec_daysanysub = prior_prec_daysanysub,
    prior_prec_pcs = prior_prec_pcs,
    prior_prec_mcs = prior_prec_mcs,
    prior_prec_cesd = prior_prec_cesd,
    prior_prec_sexrisk = prior_prec_sexrisk,
    prior_prec_drugrisk = prior_prec_drugrisk,
    prior_prec_dayslink = prior_prec_dayslink,

    sigma_upper_daysdrink = 5 * s["daysdrink"],
    sigma_upper_daysanysub = 5 * s["daysanysub"],
    sigma_upper_pcs = 5 * s["pcs"],
    sigma_upper_mcs = 5 * s["mcs"],
    sigma_upper_cesd = 5 * s["cesd"],
    sigma_upper_sexrisk = 5 * s["sexrisk"],
    sigma_upper_drugrisk = 5 * s["drugrisk"],
    sigma_upper_dayslink = 5 * s["dayslink"]
  )

  params <- c(
    "theta0",
    "beta_drink", "beta_sub", "beta_pcs", "beta_mcs",
    "beta_cesd", "beta_sexrisk", "beta_drugrisk",
    "theta90", "theta91", "theta92", "theta93", "theta94",
    "theta95", "theta96", "theta97", "theta98", "theta99",
    "sigma_daysdrink", "sigma_daysanysub", "sigma_pcs", "sigma_mcs",
    "sigma_cesd", "sigma_sexrisk", "sigma_drugrisk", "sigma_dayslink"
  )

  set.seed(seed)

  BG <- jags(
    data = datalist,
    inits = NULL,
    parameters.to.save = params,
    model.file = model_file,
    n.chains = n.chains,
    n.burnin = n.burnin,
    n.iter = n.iter,
    n.thin = n.thin,
    DIC = FALSE
  )

  BG$HELP_scaling <- scaling
  BG$HELP_source_n <- N
  BG$HELP_model_file <- model_file
  BG$HELP_treatment_probability <- mean(data$treat)
  BG
}

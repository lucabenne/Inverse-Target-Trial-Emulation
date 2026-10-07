library(readr)
library(R2jags)
library(truncnorm)

source("02_fit_HELP_bayesian_DGP.R")
source("03a_generate_HELP_synthetic_trial.R")

# -----------------------------------------------------------------------------
# 1. Source HELP data
# -----------------------------------------------------------------------------

help_raw <- read_csv("../HELP_analysis_data.csv", show_col_types = FALSE)

required <- c(
  "treat", "female", "dayslink", "daysdrink", "daysanysub",
  "pcs", "mcs", "cesd", "sexrisk", "drugrisk"
)

missing_vars <- setdiff(required, names(help_raw))
if (length(missing_vars) > 0L) {
  stop("Missing required HELP variables: ", paste(missing_vars, collapse = ", "))
}

data <- as.data.frame(help_raw[, required])

if (anyNA(data)) {
  stop("HELP validation requires complete cases for the selected variables.")
}

covariates <- c(
  "female", "daysdrink", "daysanysub", "pcs", "mcs",
  "cesd", "sexrisk", "drugrisk"
)

# -----------------------------------------------------------------------------
# 2. Fit conditional Bayesian DGP
# -----------------------------------------------------------------------------

iseed <- 101010L
set.seed(iseed)

if (!exists("BG", inherits = FALSE) || is.null(BG$HELP_source_n) || BG$HELP_source_n != nrow(data)) {
  BG <- jags_data(
    data = data,
    model_file = "../src/03_HELP_DGP_model.txt",
    seed = iseed
  )
} else {
  cat("\nUsing pre-fitted HELP Bayesian DGP from this session for validation.\n")
}


# -----------------------------------------------------------------------------
# 3. Generate one central-fit synthetic RCT
#
# Posterior means are used here to assess the central fitted DGP. In the actual
# ITTE simulation, use posterior_mode = "dataset" so each replicate receives one
# fresh joint posterior parameter draw.
# -----------------------------------------------------------------------------

set.seed(20261005L)

synthetic_rct <- generate_HELP_synthetic_trial(
  data = data,
  BG = BG,
  nsim = 1000L,
  treat_prob = mean(data$treat),
  posterior_mode = "mean"
)

cat("\nSOURCE RCT n =", nrow(data), "\n")
cat("SYNTHETIC RCT n =", nrow(synthetic_rct), "\n")
cat("Source treated proportion =", mean(data$treat), "\n")
cat("Synthetic treated proportion =", mean(synthetic_rct$treat), "\n\n")

# -----------------------------------------------------------------------------
# 4. Marginal covariate distributions
# -----------------------------------------------------------------------------

summarise_variable <- function(x) {
  c(
    mean = mean(x),
    sd = sd(x),
    median = median(x),
    q025 = as.numeric(quantile(x, 0.025)),
    q25 = as.numeric(quantile(x, 0.25)),
    q75 = as.numeric(quantile(x, 0.75)),
    q975 = as.numeric(quantile(x, 0.975))
  )
}

marginal_table <- do.call(rbind, lapply(covariates, function(v) {
  src <- summarise_variable(data[[v]])
  syn <- summarise_variable(synthetic_rct[[v]])

  pooled_sd <- sqrt((src["sd"]^2 + syn["sd"]^2) / 2)
  smd <- if (pooled_sd > 0) abs(src["mean"] - syn["mean"]) / pooled_sd else 0

  data.frame(
    variable = v,
    source_mean = src["mean"],
    synthetic_mean = syn["mean"],
    source_sd = src["sd"],
    synthetic_sd = syn["sd"],
    sd_ratio_synthetic_source = syn["sd"] / src["sd"],
    abs_standardised_difference = smd,
    source_median = src["median"],
    synthetic_median = syn["median"],
    source_q025 = src["q025"],
    synthetic_q025 = syn["q025"],
    source_q25 = src["q25"],
    synthetic_q25 = syn["q25"],
    source_q75 = src["q75"],
    synthetic_q75 = syn["q75"],
    source_q975 = src["q975"],
    synthetic_q975 = syn["q975"]
  )
}))

rownames(marginal_table) <- NULL

cat("\n================ COVARIATE MARGINALS ================\n")
print(marginal_table, row.names = FALSE, digits = 5)

cat("\nMean |SMD| =", mean(marginal_table$abs_standardised_difference), "\n")
cat("Max  |SMD| =", max(marginal_table$abs_standardised_difference), "\n")

# -----------------------------------------------------------------------------
# 5. Support checks for discrete/bounded variables
# -----------------------------------------------------------------------------

support_checks <- data.frame(
  variable = c("female", "cesd", "sexrisk", "drugrisk", "pcs", "mcs"),
  source_min = c(min(data$female), min(data$cesd), min(data$sexrisk),
                 min(data$drugrisk), min(data$pcs), min(data$mcs)),
  source_max = c(max(data$female), max(data$cesd), max(data$sexrisk),
                 max(data$drugrisk), max(data$pcs), max(data$mcs)),
  synthetic_min = c(min(synthetic_rct$female), min(synthetic_rct$cesd),
                    min(synthetic_rct$sexrisk), min(synthetic_rct$drugrisk),
                    min(synthetic_rct$pcs), min(synthetic_rct$mcs)),
  synthetic_max = c(max(synthetic_rct$female), max(synthetic_rct$cesd),
                    max(synthetic_rct$sexrisk), max(synthetic_rct$drugrisk),
                    max(synthetic_rct$pcs), max(synthetic_rct$mcs))
)

cat("\n================ SUPPORT CHECKS ================\n")
print(support_checks, row.names = FALSE)

if (any(!synthetic_rct$female %in% c(0, 1))) stop("Invalid synthetic female values.")
if (any(synthetic_rct$cesd < 0 | synthetic_rct$cesd > 60)) stop("Invalid synthetic cesd values.")
if (any(synthetic_rct$sexrisk < 0 | synthetic_rct$sexrisk > 21)) stop("Invalid synthetic sexrisk values.")
if (any(synthetic_rct$drugrisk < 0 | synthetic_rct$drugrisk > 21)) stop("Invalid synthetic drugrisk values.")

# -----------------------------------------------------------------------------
# 6. Dependence structure
# -----------------------------------------------------------------------------

source_cor <- cor(data[, covariates], use = "pairwise.complete.obs")
synthetic_cor <- cor(synthetic_rct[, covariates], use = "pairwise.complete.obs")
cor_diff <- synthetic_cor - source_cor

upper <- upper.tri(cor_diff)

correlation_table <- data.frame(
  var1 = rownames(cor_diff)[row(cor_diff)[upper]],
  var2 = colnames(cor_diff)[col(cor_diff)[upper]],
  source = source_cor[upper],
  synthetic = synthetic_cor[upper],
  difference = cor_diff[upper]
)
correlation_table$abs_difference <- abs(correlation_table$difference)
correlation_table <- correlation_table[order(-correlation_table$abs_difference), ]
rownames(correlation_table) <- NULL

cat("\n================ CORRELATIONS ================\n")
cat("Mean |correlation difference| =", mean(correlation_table$abs_difference), "\n")
cat("Max  |correlation difference| =", max(correlation_table$abs_difference), "\n\n")
print(head(correlation_table, 10), row.names = FALSE, digits = 5)

# -----------------------------------------------------------------------------
# 7. Outcome distribution
# -----------------------------------------------------------------------------

outcome_source <- summarise_variable(data$dayslink)
outcome_synthetic <- summarise_variable(synthetic_rct$dayslink)

outcome_table <- data.frame(
  dataset = c("Source RCT", "Synthetic RCT"),
  mean = c(outcome_source["mean"], outcome_synthetic["mean"]),
  sd = c(outcome_source["sd"], outcome_synthetic["sd"])
)

cat("\n================ OUTCOME DISTRIBUTION ================\n")
print(outcome_table, row.names = FALSE, digits = 5)

# -----------------------------------------------------------------------------
# 8. Treatment effect: apply the same empirical estimators to source/synthetic
# -----------------------------------------------------------------------------

raw_effect <- function(d) {
  mean(d$dayslink[d$treat == 0]) - mean(d$dayslink[d$treat == 1])
}

adjusted_effect <- function(d) {
  fit <- lm(
    dayslink ~ treat + female + daysdrink + daysanysub + pcs + mcs +
      cesd + sexrisk + drugrisk,
    data = d
  )
  -unname(coef(fit)["treat"])
}

treatment_effect_table <- data.frame(
  dataset = c("Source RCT", "Synthetic RCT"),
  raw_mean_difference_control_minus_treated = c(
    raw_effect(data),
    raw_effect(synthetic_rct)
  ),
  adjusted_linear_effect_control_minus_treated = c(
    adjusted_effect(data),
    adjusted_effect(synthetic_rct)
  )
)

cat("\n================ TREATMENT EFFECT ================\n")
print(treatment_effect_table, row.names = FALSE, digits = 5)
# theta91 is a latent-location coefficient because dayslink is lower-truncated
# at zero.  Compute the treatment contrast on the observed outcome scale.
truncnorm_lower0_mean <- function(mu, sigma) {
  z <- mu / sigma
  mills <- exp(dnorm(z, log = TRUE) - pnorm(z, log.p = TRUE))
  mu + sigma * mills
}

dgp_observed_scale_effect <- function(d, BG) {
  mode <- "mean"
  draw <- NA_integer_

  z_female <- HELP_standardize(d$female, "female", BG)
  z_daysdrink <- HELP_standardize(d$daysdrink, "daysdrink", BG)
  z_daysanysub <- HELP_standardize(d$daysanysub, "daysanysub", BG)
  z_pcs <- HELP_standardize(d$pcs, "pcs", BG)
  z_mcs <- HELP_standardize(d$mcs, "mcs", BG)
  z_cesd <- HELP_standardize(d$cesd, "cesd", BG)
  z_sexrisk <- HELP_standardize(d$sexrisk, "sexrisk", BG)
  z_drugrisk <- HELP_standardize(d$drugrisk, "drugrisk", BG)

  base_mu <-
    HELP_get_scalar(BG, "theta90", mode, draw) +
    HELP_get_scalar(BG, "theta92", mode, draw) * z_daysdrink +
    HELP_get_scalar(BG, "theta93", mode, draw) * z_daysanysub +
    HELP_get_scalar(BG, "theta94", mode, draw) * z_female +
    HELP_get_scalar(BG, "theta95", mode, draw) * z_pcs +
    HELP_get_scalar(BG, "theta96", mode, draw) * z_mcs +
    HELP_get_scalar(BG, "theta97", mode, draw) * z_cesd +
    HELP_get_scalar(BG, "theta98", mode, draw) * z_sexrisk +
    HELP_get_scalar(BG, "theta99", mode, draw) * z_drugrisk

  mu0 <- base_mu
  mu1 <- base_mu + HELP_get_scalar(BG, "theta91", mode, draw)
  sigma <- HELP_get_scalar(BG, "sigma_dayslink", mode, draw)

  mean(truncnorm_lower0_mean(mu0, sigma) -
         truncnorm_lower0_mean(mu1, sigma))
}

dgp_effect <- dgp_observed_scale_effect(data, BG)
cat("DGP implied observed-scale effect =", dgp_effect, "\n")

# -----------------------------------------------------------------------------
# 9. Save validation outputs
# -----------------------------------------------------------------------------

dir.create("results", showWarnings = FALSE, recursive = TRUE)

write.csv(marginal_table,
          "results/HELP_validation_covariate_distributions.csv",
          row.names = FALSE)
write.csv(support_checks,
          "results/HELP_validation_support.csv",
          row.names = FALSE)
write.csv(correlation_table,
          "results/HELP_validation_correlations.csv",
          row.names = FALSE)
write.csv(outcome_table,
          "results/HELP_validation_outcome.csv",
          row.names = FALSE)
write.csv(treatment_effect_table,
          "results/HELP_validation_treatment_effect.csv",
          row.names = FALSE)

save(
  BG,
  synthetic_rct,
  marginal_table,
  support_checks,
  source_cor,
  synthetic_cor,
  correlation_table,
  outcome_table,
  treatment_effect_table,
  file = "results/HELP_SYNTHETIC_RCT_VALIDATION.RData"
)

# -----------------------------------------------------------------------------
# 10. Compact final output
# -----------------------------------------------------------------------------

cat("\n\n====================================================\n")
cat("          HELP SYNTHETIC RCT VALIDATION\n")
cat("====================================================\n")
cat("Mean |SMD| =", round(mean(marginal_table$abs_standardised_difference), 5), "\n")
cat("Max  |SMD| =", round(max(marginal_table$abs_standardised_difference), 5), "\n")
cat("Mean |correlation difference| =", round(mean(correlation_table$abs_difference), 5), "\n")
cat("Max  |correlation difference| =", round(max(correlation_table$abs_difference), 5), "\n")
cat("Source raw treatment effect =", round(raw_effect(data), 3), "\n")
cat("Synthetic raw treatment effect =", round(raw_effect(synthetic_rct), 3), "\n")
cat("Source adjusted linear effect =", round(adjusted_effect(data), 3), "\n")
cat("Synthetic adjusted linear effect =", round(adjusted_effect(synthetic_rct), 3), "\n")
cat("DGP implied observed-scale effect =", round(dgp_effect, 3), "\n")
cat("====================================================\n")

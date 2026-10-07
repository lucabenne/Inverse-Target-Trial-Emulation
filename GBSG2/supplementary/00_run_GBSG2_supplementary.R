###############################################################################
# Run GBSG2 supplementary analyses
###############################################################################
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  this_file <- rstudioapi::getActiveDocumentContext()$path
  if (nzchar(this_file)) setwd(dirname(this_file))
}
dir.create("results", showWarnings = FALSE, recursive = TRUE)

source("../src/00_GBSG2_helpers.R")
source("../src/03_fit_GBSG2_covariate_DGP.R")
source("../src/05_generate_GBSG2_synthetic_trial.R")
source("../src/06_generate_additional_sampling_survival.R")
source("../src/07_generate_decoupling_covariate_survival.R")
source("../src/08_measure_mahalanobis_imbalance_survival.R")
source("../src/09_pairwise_mahalanobis_survival.R")
source("../src/10_estimate_weibull_survival_model.R")

iseed <- 100001L
set.seed(iseed)
datt <- prepare_GBSG2_data()
scaling <- get_GBSG2_survival_scaling(datt)
datt_model <- apply_GBSG2_survival_scaling(datt, scaling)
formula_full <- GBSG2_full_survival_formula()

weib <- estimate_params(datt_model, formula_full, "weibull", msd = 0, seed = iseed)
BG <- generate_BG(datt)
n_BG_draws <- GBSG2_BG_n_draws(BG)

cat("\n>>> GBSG2 synthetic-RCT validation\n")
source("01_validate_GBSG2_synthetic_RCT.R")

cat("\n>>> GBSG2 severity / overlap diagnostics\n")
source("02_GBSG2_severity_overlap_diagnostics.R")

cat("\n>>> GBSG2 fixed-final-n sensitivity\n")
source("03_GBSG2_fixed_final_n_sensitivity.R")

cat("\n>>> GBSG2 q sensitivity\n")
source("04_GBSG2_q_sensitivity.R")

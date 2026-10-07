###############################################################################
# Run HELP supplementary analyses retained for the revision
###############################################################################
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  this_file <- rstudioapi::getActiveDocumentContext()$path
  if (nzchar(this_file)) setwd(dirname(this_file))
}
dir.create("results", showWarnings = FALSE, recursive = TRUE)

source("../src/02_fit_HELP_bayesian_DGP.R")
source("../src/03a_generate_HELP_synthetic_trial.R")
source("../src/04_generate_additional_sampling.R")
source("../src/05_generate_decoupling_mahalanobis.R")
source("../src/06_generate_perfect_knowledge.R")
source("../src/07_generate_decoupling_covariate.R")
source("../src/08_measure_mahalanobis_imbalance.R")

help_raw <- readr::read_csv("../HELP_analysis_data.csv", show_col_types = FALSE)
required_vars <- c("treat","female","dayslink","daysdrink","daysanysub",
                   "pcs","mcs","cesd","sexrisk","drugrisk")
data <- as.data.frame(help_raw[, required_vars])

iseed <- 101010L
set.seed(iseed)
BG <- jags_data(data, model_file = "../src/03_HELP_DGP_model.txt", seed = iseed)

cat("\n>>> HELP synthetic-RCT validation\n")
source("../src/03b_validate_HELP_synthetic_RCT.R")

cat("\n>>> HELP severity / overlap diagnostics\n")
source("02_HELP_severity_overlap_diagnostics.R")

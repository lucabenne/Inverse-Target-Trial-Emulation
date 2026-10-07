###############################################################################
# Run final HELP analysis
###############################################################################
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  this_file <- rstudioapi::getActiveDocumentContext()$path
  if (nzchar(this_file)) setwd(dirname(this_file))
}
dir.create("results", showWarnings = FALSE, recursive = TRUE)

source("../src/02_fit_HELP_bayesian_DGP.R")
source("../src/03a_generate_HELP_synthetic_trial.R")

help_raw <- readr::read_csv("../HELP_analysis_data.csv", show_col_types = FALSE)
required_vars <- c("treat","female","dayslink","daysdrink","daysanysub",
                   "pcs","mcs","cesd","sexrisk","drugrisk")
data <- as.data.frame(help_raw[, required_vars])

iseed <- 101010L
set.seed(iseed)
BG <- jags_data(data, model_file = "../src/03_HELP_DGP_model.txt", seed = iseed)

cat("\n>>> HELP synthetic-RCT validation\n")
source("../src/03b_validate_HELP_synthetic_RCT.R")

cat("\n>>> HELP main analysis\n")
source("01_run_HELP_confounding_analysis.R")

cat("\n>>> HELP MCSE and overlap diagnostics\n")
source("17_summarize_HELP_diagnostics.R")

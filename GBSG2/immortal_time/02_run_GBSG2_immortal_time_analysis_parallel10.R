###############################################################################
# GBSG2 immortal-time-bias analysis -- 10-core parallel version
#
# Scientific design is identical to 02_run_GBSG2_immortal_time_analysis.R.
# Only the 1000 Monte Carlo replications are parallelised.
#
# Each worker processes a whole replicate (all five scenarios), while each Stan
# fit itself is kept at cores=1 in 10_estimate_weibull_survival_model.R. This
# avoids nested oversubscription and targets approximately 10 concurrent cores.
###############################################################################

library(TH.data)
library(survival)
library(R2jags)
library(truncnorm)
library(parallel)

source("../src/00_GBSG2_helpers.R")
source("../src/03_fit_GBSG2_covariate_DGP.R")
source("../src/05_generate_GBSG2_synthetic_trial.R")
source("../src/10_estimate_weibull_survival_model.R")
source("../src/11_adjust_PTDM.R")

iseed <- 100001L
n_rep <- 1000L
nsim <- 1000L
requested_cores <- 10L

target_reclassified <- c(5, 10, 15, 20, 25)
immortal_time_variance <- 20
calibration_reps <- 1000L

set.seed(iseed)

datt <- prepare_GBSG2_data()
scaling <- get_GBSG2_survival_scaling(datt)
datt_model <- apply_GBSG2_survival_scaling(datt, scaling)
formula_full <- GBSG2_full_survival_formula()
distrib <- "weibull"

cat("\n>>> FIT SOURCE WEIBULL EVENT/CENSORING DGP\n")
weib <- estimate_params(
  datt_model,
  formula_full,
  distrib,
  msd = 0,
  seed = iseed
)
cat("Source adjusted Weibull MSD =", weib$mean_surv_diff, "\n")

cat("\n>>> FIT GBSG2 COVARIATE DGP AND GENERATE BASELINE SYNTHETIC RCT\n")
BG <- generate_BG(datt)
set.seed(iseed)
covariates <- generate_GBSG2_covariates(datt, BG = BG, nsim = nsim)
treat <- rbinom(nsim, 1, 0.5)
sim_base <- data.frame(treat = treat, covariates)

# Arrival-process variables retained unchanged from the previous analysis.
n_hospitals <- 5L
lambda <- numeric(n_hospitals)
for(i in seq_len(n_hospitals)){
  lambda[i] <- rgamma(1, shape = 1, scale = 1)
  while(lambda[i] < 0.01){
    lambda[i] <- rgamma(1, shape = 1, scale = 0.01)
  }
}

hospital <- sample(seq_len(n_hospitals), nrow(sim_base), replace = TRUE)
arr_time <- numeric(nrow(sim_base))
for(i in seq_len(nrow(sim_base))){
  arr_time[i] <- rexp(1, lambda[hospital[i]])
}
for(i in seq_len(n_hospitals)){
  idx <- which(hospital == i)
  arr_time[idx] <- cumsum(arr_time[idx])
}
arr_time <- round(arr_time)

idx_assigned_treated <- which(sim_base$treat == 1)
treated_profiles <- sim_base[idx_assigned_treated, , drop = FALSE]

###############################################################################
# Automatic calibration to approximately 5, 10, 15, 20, 25 reclassified.
###############################################################################

cat("\n>>> CALIBRATING IMMORTAL-TIME SCENARIOS\n")
set.seed(iseed + 7000L)

observed0_calibration <- matrix(
  NA_real_,
  nrow = calibration_reps,
  ncol = length(idx_assigned_treated)
)

for(b in seq_len(calibration_reps)){
  event0_b <- GBSG2_draw_weibull_time(
    treated_profiles,
    weib$beta_event,
    weib$alpha_event,
    scaling,
    treat_override = 0
  )
  cens0_b <- GBSG2_draw_weibull_time(
    treated_profiles,
    weib$beta_cens,
    weib$alpha_cens,
    scaling,
    treat_override = 0
  )
  observed0_calibration[b, ] <- pmin(event0_b, cens0_b)
}

expected_reclassified <- function(mean_immortal_time){
  if(!is.finite(mean_immortal_time) || mean_immortal_time <= 0){
    return(0)
  }
  shape <- mean_immortal_time^2 / immortal_time_variance
  rate <- mean_immortal_time / immortal_time_variance
  p_reclassified <- pgamma(
    observed0_calibration,
    shape = shape,
    rate = rate,
    lower.tail = FALSE
  )
  mean(rowSums(p_reclassified))
}

calibrate_mean <- function(target){
  f <- function(m) expected_reclassified(m) - target
  lo <- 0.10
  hi <- 50
  while(f(hi) < 0){
    hi <- hi * 1.5
    if(hi > 5000) stop("Unable to bracket calibration target: ", target)
  }
  uniroot(f, interval = c(lo, hi), tol = 1e-6)$root
}

mean_immortal_time <- vapply(target_reclassified, calibrate_mean, numeric(1))
immortal_shape <- mean_immortal_time^2 / immortal_time_variance
immortal_rate <- mean_immortal_time / immortal_time_variance
calibrated_expected_reclassified <- vapply(
  mean_immortal_time,
  expected_reclassified,
  numeric(1)
)

calibration_table <- data.frame(
  scenario = seq_along(target_reclassified),
  target_reclassified = target_reclassified,
  calibrated_mean_immortal_time = mean_immortal_time,
  gamma_shape = immortal_shape,
  gamma_rate = immortal_rate,
  gamma_variance = immortal_time_variance,
  expected_reclassified = calibrated_expected_reclassified
)

cat("\n================ CALIBRATED SCENARIOS ================\n")
print(calibration_table, row.names = FALSE, digits = 6)
write.csv(
  calibration_table,
  "GBSG2_immortal_time_calibration.csv",
  row.names = FALSE
)

K <- length(target_reclassified)

run_one_immortal_rep <- function(h){
  set.seed(iseed + h)

  event0 <- GBSG2_draw_weibull_time(
    sim_base,
    weib$beta_event,
    weib$alpha_event,
    scaling,
    treat_override = 0
  )
  cens0 <- GBSG2_draw_weibull_time(
    sim_base,
    weib$beta_cens,
    weib$alpha_cens,
    scaling,
    treat_override = 0
  )

  observed0 <- pmin(event0, cens0)
  status0 <- as.numeric(event0 < cens0)

  event1 <- GBSG2_draw_weibull_time(
    treated_profiles,
    weib$beta_event,
    weib$alpha_event,
    scaling,
    treat_override = 1
  )
  cens1 <- GBSG2_draw_weibull_time(
    treated_profiles,
    weib$beta_cens,
    weib$alpha_cens,
    scaling,
    treat_override = 1
  )

  out_naive <- rep(NA_real_, K)
  out_pstdm <- rep(NA_real_, K)
  out_reclassified <- rep(NA_real_, K)

  for(k in seq_len(K)){
    imm_treated <- rgamma(
      length(idx_assigned_treated),
      shape = immortal_shape[k],
      rate = immortal_rate[k]
    )

    immortal_time <- rep(0, nrow(sim_base))
    immortal_time[idx_assigned_treated] <- imm_treated

    time <- observed0
    event <- status0
    new_treat <- rep(0, nrow(sim_base))

    eligible_local <- which(
      imm_treated < observed0[idx_assigned_treated]
    )

    if(length(eligible_local)){
      eligible_global <- idx_assigned_treated[eligible_local]
      post_obs <- pmin(event1[eligible_local], cens1[eligible_local])
      time[eligible_global] <- imm_treated[eligible_local] + post_obs
      event[eligible_global] <- as.numeric(
        event1[eligible_local] < cens1[eligible_local]
      )
      new_treat[eligible_global] <- 1
    }

    out_reclassified[k] <- length(idx_assigned_treated) - length(eligible_local)

    sim_data <- cbind(
      sim_base,
      hospital = hospital,
      arr_time = arr_time,
      time = time,
      event = event,
      immortal_time = immortal_time,
      new_treat = new_treat
    )

    sim_data1 <- pstdm(sim_data)
    sim_data$treat <- sim_data$new_treat

    model_data <- apply_GBSG2_survival_scaling(sim_data, scaling)
    model_data1 <- apply_GBSG2_survival_scaling(sim_data1, scaling)

    fit_naive <- estimate_params(
      model_data,
      formula_full,
      distrib,
      msd = 1,
      seed = iseed + h * 100L + k * 2L
    )

    fit_pstdm <- estimate_params(
      model_data1,
      formula_full,
      distrib,
      msd = 1,
      seed = iseed + h * 100L + k * 2L + 1L
    )

    out_naive[k] <- fit_naive$mean_surv_diff
    out_pstdm[k] <- fit_pstdm$mean_surv_diff
  }

  list(
    replicate = h,
    naive = out_naive,
    ptdm = out_pstdm,
    reclassified = out_reclassified
  )
}

available_cores <- parallel::detectCores(logical = TRUE)
if(is.na(available_cores)) available_cores <- requested_cores
n_cores <- max(1L, min(requested_cores, available_cores))

cat("\n>>> PARALLEL IMMORTAL-TIME MONTE CARLO\n")
cat("Replicates:", n_rep, "\n")
cat("Workers:", n_cores, "\n")
cat("Each worker runs complete replicates; Stan itself remains cores=1 per fit.\n\n")

if(.Platform$OS.type != "windows"){
  # Forking is efficient on macOS/Linux because the large fitted objects can be
  # inherited copy-on-write by workers.
  results <- parallel::mclapply(
    seq_len(n_rep),
    run_one_immortal_rep,
    mc.cores = n_cores,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
} else {
  # Cross-platform fallback for Windows.
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterEvalQ(cl, {
    library(TH.data)
    library(survival)
    library(R2jags)
    library(truncnorm)
    source("../src/00_GBSG2_helpers.R")
    source("../src/10_estimate_weibull_survival_model.R")
    source("../src/11_adjust_PTDM.R")
    NULL
  })

  parallel::clusterExport(
    cl,
    varlist = c(
      "iseed", "K", "sim_base", "weib", "scaling", "formula_full", "distrib",
      "treated_profiles", "idx_assigned_treated", "hospital", "arr_time",
      "immortal_shape", "immortal_rate", "run_one_immortal_rep"
    ),
    envir = environment()
  )

  results <- parallel::parLapply(
    cl,
    seq_len(n_rep),
    run_one_immortal_rep
  )
}

# Reassemble in replicate order regardless of completion order.
ord <- order(vapply(results, `[[`, integer(1), "replicate"))
results <- results[ord]

mean_sim <- as.data.frame(do.call(rbind, lapply(results, `[[`, "naive")))
mean_pstdm <- as.data.frame(do.call(rbind, lapply(results, `[[`, "ptdm")))
reclassified_ind <- as.data.frame(do.call(rbind, lapply(results, `[[`, "reclassified")))

colnames(mean_sim) <- colnames(mean_pstdm) <- colnames(reclassified_ind) <-
  paste0("scenario", seq_len(K))

# Backward-compatible alias for old plotting code.
sw_ind <- reclassified_ind

immortal_summary <- data.frame(
  scenario = seq_len(K),
  target_reclassified_individuals = target_reclassified,
  shape = immortal_shape,
  rate = immortal_rate,
  mean_immortal_time = mean_immortal_time,
  source_adjusted_MSD = rep(weib$mean_surv_diff, K),
  mean_simulated_MSD = colMeans(mean_sim, na.rm = TRUE),
  mean_PTDM_MSD = colMeans(mean_pstdm, na.rm = TRUE),
  mean_reclassified_individuals = colMeans(reclassified_ind, na.rm = TRUE),
  mean_switched_individuals = colMeans(reclassified_ind, na.rm = TRUE)
)

cat("\n================ IMMORTAL TIME SUMMARY ================\n")
print(immortal_summary, row.names = FALSE, digits = 5)

write.csv(
  immortal_summary,
  "GBSG2_immortal_time_summary_weibull.csv",
  row.names = FALSE
)

save(
  weib,
  BG,
  target_reclassified,
  immortal_time_variance,
  immortal_shape,
  immortal_rate,
  mean_immortal_time,
  calibration_table,
  mean_sim,
  mean_pstdm,
  reclassified_ind,
  sw_ind,
  immortal_summary,
  file = "GBSG2_immortal_time_results_weibull.RData"
)

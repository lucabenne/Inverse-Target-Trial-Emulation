###############################################################################
# GBSG2 immortal-time-bias analysis -- canonical sequential version
#
# Five scenarios are calibrated to yield approximately 5, 10, 15, 20 and 25
# individuals reclassified as untreated because treatment initiation occurs
# after their untreated event/censoring time.
#
# IMPORTANT:
# - The reclassified count is an interpretable proxy for scenario severity.
# - Immortal-time bias is generated more generally by assigning pre-treatment
#   follow-up to subjects subsequently classified as treated; it is not confined
#   to the reclassified individuals themselves.
# - The Gamma treatment-initiation-time distribution keeps variance fixed at 20
#   while its mean is calibrated automatically for each target count.
###############################################################################

library(TH.data)
library(survival)
library(R2jags)
library(truncnorm)

source("../src/00_GBSG2_helpers.R")
source("../src/03_fit_GBSG2_covariate_DGP.R")
source("../src/05_generate_GBSG2_synthetic_trial.R")
source("../src/10_estimate_weibull_survival_model.R")
source("../src/11_adjust_PTDM.R")

iseed <- 100001L
n_rep <- 1000L
nsim <- 1000L

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

# Baseline synthetic RCT: identical architecture to the validated pipeline.
cat("\n>>> FIT GBSG2 COVARIATE DGP AND GENERATE BASELINE SYNTHETIC RCT\n")
BG <- generate_BG(datt)
set.seed(iseed)
covariates <- generate_GBSG2_covariates(datt, BG = BG, nsim = nsim)
treat <- rbinom(nsim, 1, 0.5)
sim_base <- data.frame(treat = treat, covariates)

# Arrival-process variables retained unchanged from the previous immortal-time run.
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
# Calibrate treatment-initiation-time means to target reclassified counts.
#
# For Gamma(shape, rate): mean = shape/rate and variance = shape/rate^2.
# Holding variance V fixed gives:
#   shape = mean^2 / V
#   rate  = mean / V
#
# For a subject with untreated observed time t0,
# P(reclassified) = P(initiation_time >= t0).
# The calibration averages this probability over repeated untreated event/
# censoring draws from the fitted Weibull DGP. No survival-model refits occur.
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

  # Probability that treatment initiation happens after untreated event/censoring.
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

###############################################################################
# Main Monte Carlo analysis
###############################################################################

K <- length(target_reclassified)
mean_sim <- matrix(NA_real_, nrow = n_rep, ncol = K)
mean_pstdm <- matrix(NA_real_, nrow = n_rep, ncol = K)
reclassified_ind <- matrix(NA_real_, nrow = n_rep, ncol = K)

colnames(mean_sim) <- colnames(mean_pstdm) <- colnames(reclassified_ind) <-
  paste0("scenario", seq_len(K))

for(h in seq_len(n_rep)){
  set.seed(iseed + h)

  if(h == 1L || h %% 10L == 0L){
    cat("Immortal-time replicate", h, "of", n_rep, "\n")
  }

  # Before actual treatment initiation, everyone follows the untreated DGP.
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

  # Common post-treatment draws across all five scenarios in the replicate.
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

    reclassified <- length(idx_assigned_treated) - length(eligible_local)
    reclassified_ind[h, k] <- reclassified

    sim_data <- cbind(
      sim_base,
      hospital = hospital,
      arr_time = arr_time,
      time = time,
      event = event,
      immortal_time = immortal_time,
      new_treat = new_treat
    )

    # PTDM correction retained unchanged.
    sim_data1 <- pstdm(sim_data)

    # Naive analysis treats eventual treatment status as a baseline exposure.
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

    mean_sim[h, k] <- fit_naive$mean_surv_diff
    mean_pstdm[h, k] <- fit_pstdm$mean_surv_diff
  }
}

mean_sim <- as.data.frame(mean_sim)
mean_pstdm <- as.data.frame(mean_pstdm)
reclassified_ind <- as.data.frame(reclassified_ind)

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
  # Legacy column name retained so existing plotting scripts do not break.
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

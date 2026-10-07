###############################################################################
# HELP Monte Carlo errors and confounding/overlap diagnostics
# Requires 01_run_HELP_confounding_analysis.R to have been run in this session.
###############################################################################

covariate_names <- c(
  "female", "daysdrink", "daysanysub", "pcs", "mcs",
  "cesd", "sexrisk", "drugrisk"
)

mc_summary <- function(x, target) {
  if (length(target) == 1L) target <- rep(target, length(x))
  if (length(target) != length(x)) stop("Target and estimate vectors must have the same length.")
  ok <- is.finite(x) & is.finite(target)
  err <- abs(x[ok] - target[ok])
  c(
    mean_target = mean(target[ok]),
    sd_target = sd(target[ok]),
    mean_estimate = mean(x[ok]),
    empirical_sd = sd(x[ok]),
    mcse_mean = sd(x[ok]) / sqrt(sum(ok)),
    mae = mean(err),
    mcse_mae = sd(err) / sqrt(length(err))
  )
}

scenario_vectors <- list(
  "Additional Sampling" = list(
    Unadjusted = ate_before2,
    Mahalanobis = att_tte2,
    IPW = ate_ipw2,
    Regression = ate_reg2,
    AIPW = ate_aipw2
  ),
  "Decoupling Mahalanobis" = list(
    Unadjusted = ate_before3,
    Mahalanobis = att_tte3,
    IPW = ate_ipw3,
    Regression = ate_reg3,
    AIPW = ate_aipw3
  ),
  "Perfect Knowledge" = list(
    Unadjusted = ate_before4,
    Mahalanobis = att_tte4,
    IPW = ate_ipw4,
    Regression = ate_reg4,
    AIPW = ate_aipw4
  ),
  "Decoupling Covariate" = list(
    Unadjusted = ate_before5,
    Mahalanobis = att_tte5,
    IPW = ate_ipw5,
    Regression = ate_reg5,
    AIPW = ate_aipw5
  )
)

scenario_targets <- list(
  "Additional Sampling" = target2,
  "Decoupling Mahalanobis" = target3,
  "Perfect Knowledge" = target4,
  "Decoupling Covariate" = target5
)

mc_rows <- list()
k <- 1L
for (sc in names(scenario_vectors)) {
  for (est in names(scenario_vectors[[sc]])) {
    sm <- mc_summary(scenario_vectors[[sc]][[est]], scenario_targets[[sc]])
    mc_rows[[k]] <- data.frame(
      scenario = sc,
      estimator = est,
      mean_target = sm["mean_target"],
      sd_target = sm["sd_target"],
      mean_estimate = sm["mean_estimate"],
      empirical_sd = sm["empirical_sd"],
      mcse_mean = sm["mcse_mean"],
      mae = sm["mae"],
      mcse_mae = sm["mcse_mae"],
      row.names = NULL
    )
    k <- k + 1L
  }
}
mc_table2 <- do.call(rbind, mc_rows)
write.csv(mc_table2, "results/HELP_table2_MCSE.csv", row.names = FALSE)

misspec_targets <- scenario_targets

misspec_vectors <- list(
  "Outcome misspecification" = list(
    "Additional Sampling" = list(RB = ate_reg2_biased_outcome, AIPW = ate_aipw2_biased_outcome),
    "Decoupling Mahalanobis" = list(RB = ate_reg3_biased_outcome, AIPW = ate_aipw3_biased_outcome),
    "Perfect Knowledge" = list(RB = ate_reg4_biased_outcome, AIPW = ate_aipw4_biased_outcome),
    "Decoupling Covariate" = list(RB = ate_reg5_biased_outcome, AIPW = ate_aipw5_biased_outcome)
  ),
  "Propensity-score misspecification" = list(
    "Additional Sampling" = list(IPW = ate_ipw2_biased_prop, AIPW = ate_aipw2_biased_prop),
    "Decoupling Mahalanobis" = list(IPW = ate_ipw3_biased_prop, AIPW = ate_aipw3_biased_prop),
    "Perfect Knowledge" = list(IPW = ate_ipw4_biased_prop, AIPW = ate_aipw4_biased_prop),
    "Decoupling Covariate" = list(IPW = ate_ipw5_biased_prop, AIPW = ate_aipw5_biased_prop)
  )
)

mc3_rows <- list()
k <- 1L
for (analysis_name in names(misspec_vectors)) {
  for (sc in names(misspec_vectors[[analysis_name]])) {
    for (est in names(misspec_vectors[[analysis_name]][[sc]])) {
      x <- misspec_vectors[[analysis_name]][[sc]][[est]]
      sm <- mc_summary(x, misspec_targets[[sc]])
      mc3_rows[[k]] <- data.frame(
        analysis = analysis_name,
        scenario = sc,
        estimator = est,
        mean_target = sm["mean_target"],
        sd_target = sm["sd_target"],
        mean_estimate = sm["mean_estimate"],
        empirical_sd = sm["empirical_sd"],
        mcse_mean = sm["mcse_mean"],
        mae = sm["mae"],
        mcse_mae = sm["mcse_mae"],
        row.names = NULL
      )
      k <- k + 1L
    }
  }
}
mc_table3 <- do.call(rbind, mc3_rows)
write.csv(mc_table3, "results/HELP_table3_MCSE.csv", row.names = FALSE)

ps_overlap_coefficient <- function(ps, z) {
  ps1 <- ps[z == 1]
  ps0 <- ps[z == 0]
  if (length(ps1) < 2L || length(ps0) < 2L) return(NA_real_)
  d1 <- density(ps1, from = 0, to = 1, n = 512)
  d0 <- density(ps0, from = 0, to = 1, n = 512)
  dx <- d1$x[2] - d1$x[1]
  sum(pmin(d1$y, d0$y)) * dx
}

diagnose_dataset <- function(d) {
  z <- d$treat
  x <- d[, covariate_names, drop = FALSE]

  smd <- vapply(covariate_names, function(v) {
    x1 <- d[[v]][z == 1]
    x0 <- d[[v]][z == 0]
    denom <- sqrt((var(x1) + var(x0)) / 2)
    if (!is.finite(denom) || denom == 0) return(0)
    (mean(x1) - mean(x0)) / denom
  }, numeric(1))

  ps_fit <- suppressWarnings(glm(
    treat ~ female + daysdrink + daysanysub + pcs + mcs + cesd + sexrisk + drugrisk,
    family = binomial,
    data = d
  ))
  ps <- as.numeric(predict(ps_fit, type = "response"))
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  w <- ifelse(z == 1, 1 / ps, 1 / (1 - ps))
  ess <- (sum(w)^2) / sum(w^2)

  c(
    n = nrow(d),
    treated_proportion = mean(z),
    mahalanobis = tryCatch(overlapping(d, 3), error = function(e) NA_real_),
    mean_abs_smd = mean(abs(smd)),
    max_abs_smd = max(abs(smd)),
    ps_overlap = ps_overlap_coefficient(ps, z),
    ess_ratio = ess / nrow(d),
    max_weight = max(w)
  )
}

scenario_data <- list(
  "Additional Sampling" = list_data2,
  "Decoupling Mahalanobis" = list_data3,
  "Perfect Knowledge" = list_data4,
  "Decoupling Covariate" = list_data5
)

severity_rep <- list()
severity_summary_rows <- list()
k <- 1L
for (sc in names(scenario_data)) {
  cat("Severity diagnostics:", sc, "\n")
  m <- t(vapply(scenario_data[[sc]], diagnose_dataset, numeric(8)))
  m <- as.data.frame(m)
  m$replicate <- seq_len(nrow(m))
  m$scenario <- sc
  severity_rep[[sc]] <- m

  severity_summary_rows[[k]] <- data.frame(
    scenario = sc,
    mean_n = mean(m$n, na.rm = TRUE),
    mean_treated_proportion = mean(m$treated_proportion, na.rm = TRUE),
    mean_mahalanobis = mean(m$mahalanobis, na.rm = TRUE),
    mean_abs_smd = mean(m$mean_abs_smd, na.rm = TRUE),
    mean_max_abs_smd = mean(m$max_abs_smd, na.rm = TRUE),
    mean_ps_overlap = mean(m$ps_overlap, na.rm = TRUE),
    mean_ess_ratio = mean(m$ess_ratio, na.rm = TRUE),
    median_max_weight = median(m$max_weight, na.rm = TRUE),
    p95_max_weight = as.numeric(quantile(m$max_weight, 0.95, na.rm = TRUE)),
    row.names = NULL
  )
  k <- k + 1L
}

severity_summary <- do.call(rbind, severity_summary_rows)
severity_replicates <- do.call(rbind, severity_rep)
write.csv(severity_summary, "results/HELP_severity_overlap_summary.csv", row.names = FALSE)
write.csv(severity_replicates, "results/HELP_severity_overlap_replicates.csv", row.names = FALSE)

save(
  mc_table2, mc_table3, severity_summary, severity_replicates,
  file = "results/HELP_revision_diagnostics.RData"
)

cat("\nMonte Carlo error and severity/overlap diagnostics saved in ./results/\n")
print(severity_summary, row.names = FALSE)

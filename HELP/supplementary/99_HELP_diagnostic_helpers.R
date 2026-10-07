###############################################################################
# 99_HELP_diagnostic_helpers.R
###############################################################################

HELP_covariates_diag <- c(
  "female","daysdrink","daysanysub","pcs","mcs","cesd","sexrisk","drugrisk"
)

HELP_safe_smd <- function(x1, x0) {
  den <- sqrt((var(x1) + var(x0)) / 2)
  if (!is.finite(den) || den == 0) return(0)
  (mean(x1) - mean(x0)) / den
}

HELP_ps_overlap_coefficient <- function(ps, z) {
  ps1 <- ps[z == 1]
  ps0 <- ps[z == 0]
  if (length(ps1) < 2L || length(ps0) < 2L) return(NA_real_)
  d1 <- density(ps1, from = 0, to = 1, n = 512)
  d0 <- density(ps0, from = 0, to = 1, n = 512)
  dx <- d1$x[2] - d1$x[1]
  sum(pmin(d1$y, d0$y)) * dx
}

HELP_diag_one <- function(d, source_data = NULL, target = NA_real_) {
  z <- d$treat

  smd <- vapply(HELP_covariates_diag, function(v) {
    HELP_safe_smd(d[[v]][z == 1], d[[v]][z == 0])
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

  pooled_mean_abs_smd <- NA_real_
  pooled_max_abs_smd <- NA_real_
  if (!is.null(source_data)) {
    pooled_smd <- vapply(HELP_covariates_diag, function(v) {
      x0 <- source_data[[v]]
      x1 <- d[[v]]
      den <- sqrt((var(x0) + var(x1)) / 2)
      if (!is.finite(den) || den == 0) return(0)
      (mean(x1) - mean(x0)) / den
    }, numeric(1))
    pooled_mean_abs_smd <- mean(abs(pooled_smd))
    pooled_max_abs_smd <- max(abs(pooled_smd))
  }

  crude <- mean(d$dayslink[d$treat == 0]) - mean(d$dayslink[d$treat == 1])

  out <- c(
    n = nrow(d),
    treated_proportion = mean(z),
    mahalanobis = tryCatch(overlapping(d, 3), error = function(e) NA_real_),
    mean_abs_smd = mean(abs(smd)),
    max_abs_smd = max(abs(smd)),
    ps_overlap = HELP_ps_overlap_coefficient(ps, z),
    ess_ratio = ess / nrow(d),
    max_weight = max(w),
    pooled_mean_abs_smd = pooled_mean_abs_smd,
    pooled_max_abs_smd = pooled_max_abs_smd,
    crude_effect = crude,
    crude_abs_target_difference = if (is.finite(target)) abs(crude - target) else NA_real_
  )

  attr(out, "smd") <- smd
  out
}

HELP_diag_list <- function(L, scenario, source_data, target) {
  m <- t(vapply(L, HELP_diag_one, numeric(12),
                source_data = source_data, target = target))
  m <- as.data.frame(m)
  m$replicate <- seq_len(nrow(m))
  m$scenario <- scenario

  smd_long <- do.call(rbind, lapply(seq_along(L), function(i) {
    one <- HELP_diag_one(L[[i]], source_data = source_data, target = target)
    s <- attr(one, "smd")
    data.frame(
      scenario = scenario,
      replicate = i,
      covariate = names(s),
      smd = as.numeric(s),
      abs_smd = abs(as.numeric(s)),
      row.names = NULL
    )
  }))

  list(replicates = m, smd = smd_long)
}

HELP_summarise_diag <- function(m) {
  data.frame(
    scenario = m$scenario[1],
    n_replicates = nrow(m),
    mean_n = mean(m$n, na.rm = TRUE),
    mean_treated_proportion = mean(m$treated_proportion, na.rm = TRUE),
    mean_mahalanobis = mean(m$mahalanobis, na.rm = TRUE),
    mean_abs_smd = mean(m$mean_abs_smd, na.rm = TRUE),
    mean_max_abs_smd = mean(m$max_abs_smd, na.rm = TRUE),
    mean_ps_overlap = mean(m$ps_overlap, na.rm = TRUE),
    mean_ess_ratio = mean(m$ess_ratio, na.rm = TRUE),
    median_max_weight = median(m$max_weight, na.rm = TRUE),
    p95_max_weight = unname(quantile(m$max_weight, .95, na.rm = TRUE)),
    mean_pooled_abs_smd = mean(m$pooled_mean_abs_smd, na.rm = TRUE),
    mean_pooled_max_abs_smd = mean(m$pooled_max_abs_smd, na.rm = TRUE),
    mean_crude_effect = mean(m$crude_effect, na.rm = TRUE),
    mean_crude_abs_target_difference = mean(m$crude_abs_target_difference, na.rm = TRUE),
    row.names = NULL
  )
}

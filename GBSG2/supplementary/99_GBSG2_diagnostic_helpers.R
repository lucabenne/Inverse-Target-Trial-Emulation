###############################################################################
# 99_GBSG2_diagnostic_helpers.R
###############################################################################

GBSG2_diag_covariates <- c(
  "age","menostat","tsize","tgrade","pnodes","progrec","estrec"
)

GBSG2_safe_smd <- function(x1, x0) {
  den <- sqrt((var(x1) + var(x0)) / 2)
  if (!is.finite(den) || den == 0) return(0)
  (mean(x1) - mean(x0)) / den
}

GBSG2_ps_overlap_coefficient <- function(ps, z) {
  ps1 <- ps[z == 1]
  ps0 <- ps[z == 0]
  if (length(ps1) < 2L || length(ps0) < 2L) return(NA_real_)
  d1 <- density(ps1, from = 0, to = 1, n = 512)
  d0 <- density(ps0, from = 0, to = 1, n = 512)
  dx <- d1$x[2] - d1$x[1]
  sum(pmin(d1$y, d0$y)) * dx
}

GBSG2_mahalanobis_all_covariates <- function(d) {
  z <- d$treat
  x1 <- d[z == 1, GBSG2_diag_covariates, drop = FALSE]
  x0 <- d[z == 0, GBSG2_diag_covariates, drop = FALSE]
  S <- (cov(x1) + cov(x0)) / 2
  delta <- colMeans(x1) - colMeans(x0)
  as.numeric(sqrt(t(delta) %*% solve(S) %*% delta))
}

GBSG2_diag_one <- function(d, source_data = NULL, mahal_index = NULL) {
  z <- d$treat

  smd <- vapply(GBSG2_diag_covariates, function(v) {
    GBSG2_safe_smd(d[[v]][z == 1], d[[v]][z == 0])
  }, numeric(1))

  ps_fit <- suppressWarnings(glm(
    treat ~ age + menostat + tsize + tgrade + pnodes + progrec + estrec,
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
    pooled_smd <- vapply(GBSG2_diag_covariates, function(v) {
      x0 <- source_data[[v]]
      x1 <- d[[v]]
      den <- sqrt((var(x0) + var(x1)) / 2)
      if (!is.finite(den) || den == 0) return(0)
      (mean(x1) - mean(x0)) / den
    }, numeric(1))
    pooled_mean_abs_smd <- mean(abs(pooled_smd))
    pooled_max_abs_smd <- max(abs(pooled_smd))
  }

  mahal_original <- if (is.null(mahal_index)) {
    NA_real_
  } else {
    tryCatch(overlapping_surv(d, mahal_index), error = function(e) NA_real_)
  }

  out <- c(
    n = nrow(d),
    treated_proportion = mean(z),
    event_proportion = mean(d$event),
    mahalanobis_original = mahal_original,
    mahalanobis_all_covariates =
      tryCatch(GBSG2_mahalanobis_all_covariates(d), error = function(e) NA_real_),
    mean_abs_smd = mean(abs(smd)),
    max_abs_smd = max(abs(smd)),
    ps_overlap = GBSG2_ps_overlap_coefficient(ps, z),
    ess_ratio = ess / nrow(d),
    max_weight = max(w),
    pooled_mean_abs_smd = pooled_mean_abs_smd,
    pooled_max_abs_smd = pooled_max_abs_smd
  )
  attr(out, "smd") <- smd
  out
}

GBSG2_summarise_replicates <- function(m, extra = list()) {
  out <- data.frame(
    n_replicates = nrow(m),
    mean_n = mean(m$n, na.rm = TRUE),
    mean_treated_proportion = mean(m$treated_proportion, na.rm = TRUE),
    mean_event_proportion = mean(m$event_proportion, na.rm = TRUE),
    mean_mahalanobis_original = mean(m$mahalanobis_original, na.rm = TRUE),
    mean_mahalanobis_all_covariates = mean(m$mahalanobis_all_covariates, na.rm = TRUE),
    mean_abs_smd = mean(m$mean_abs_smd, na.rm = TRUE),
    mean_max_abs_smd = mean(m$max_abs_smd, na.rm = TRUE),
    mean_ps_overlap = mean(m$ps_overlap, na.rm = TRUE),
    mean_ess_ratio = mean(m$ess_ratio, na.rm = TRUE),
    median_max_weight = median(m$max_weight, na.rm = TRUE),
    p95_max_weight = unname(quantile(m$max_weight, .95, na.rm = TRUE)),
    mean_pooled_abs_smd = mean(m$pooled_mean_abs_smd, na.rm = TRUE),
    mean_pooled_max_abs_smd = mean(m$pooled_max_abs_smd, na.rm = TRUE),
    row.names = NULL
  )
  if (length(extra)) {
    for (nm in names(extra)) out[[nm]] <- extra[[nm]]
  }
  out
}

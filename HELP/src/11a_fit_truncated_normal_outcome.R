# Correctly specified lower-truncated Normal outcome nuisance model for HELP.
# The DGP is Y | A,X ~ Normal(mu, sigma^2) truncated below at 0.
# This helper estimates the latent-location regression by maximum likelihood and
# returns predictions on the observed outcome scale E[Y | A,X,Y>0].

HELP_fit_truncnorm_outcome <- function(data) {
  covars <- c("female","daysdrink","daysanysub","pcs","mcs","cesd","sexrisk","drugrisk")
  form <- as.formula(paste("dayslink ~ treat +", paste(covars, collapse = " + ")))

  X_raw <- model.matrix(form, data = data)
  y <- as.numeric(data$dayslink)
  if (any(!is.finite(y)) || any(y <= 0)) {
    stop("Truncated-Normal outcome fit requires strictly positive finite dayslink values.")
  }

  # Standardize non-intercept design columns internally for numerical stability.
  x_center <- rep(0, ncol(X_raw))
  x_scale <- rep(1, ncol(X_raw))
  names(x_center) <- names(x_scale) <- colnames(X_raw)
  non_intercept <- which(colnames(X_raw) != "(Intercept)")
  if (length(non_intercept) > 0L) {
    x_center[non_intercept] <- colMeans(X_raw[, non_intercept, drop = FALSE])
    x_scale[non_intercept] <- apply(X_raw[, non_intercept, drop = FALSE], 2L, sd)
    x_scale[!is.finite(x_scale) | x_scale == 0] <- 1
  }
  X <- sweep(X_raw, 2L, x_center, "-")
  X <- sweep(X, 2L, x_scale, "/")

  lm_start <- lm.fit(X, y)
  beta0 <- lm_start$coefficients
  beta0[!is.finite(beta0)] <- 0
  resid0 <- as.numeric(y - X %*% beta0)
  sigma0 <- sd(resid0)
  if (!is.finite(sigma0) || sigma0 <= 0) sigma0 <- sd(y)
  if (!is.finite(sigma0) || sigma0 <= 0) sigma0 <- 1

  nll <- function(par) {
    beta <- par[seq_len(ncol(X))]
    sigma <- exp(par[ncol(X) + 1L])
    mu <- as.numeric(X %*% beta)
    z <- mu / sigma

    ll <- dnorm(y, mean = mu, sd = sigma, log = TRUE) -
      pnorm(z, log.p = TRUE)

    if (any(!is.finite(ll))) return(.Machine$double.xmax / 100)
    -sum(ll)
  }

  start <- c(beta0, log(sigma0))
  fit <- optim(start, nll, method = "BFGS", control = list(maxit = 1000, reltol = 1e-9))

  if (fit$convergence != 0L || !is.finite(fit$value)) {
    fit <- optim(start, nll, method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-9))
  }

  if (fit$convergence != 0L || !is.finite(fit$value)) {
    stop("Lower-truncated Normal outcome regression did not converge.")
  }

  beta <- fit$par[seq_len(ncol(X))]
  names(beta) <- colnames(X)
  sigma <- exp(fit$par[ncol(X) + 1L])

  structure(
    list(
      beta = beta,
      sigma = sigma,
      formula = form,
      x_center = x_center,
      x_scale = x_scale,
      convergence = fit$convergence,
      value = fit$value
    ),
    class = "HELP_truncnorm_fit"
  )
}

HELP_predict_truncnorm_mean <- function(fit, newdata) {
  if (!inherits(fit, "HELP_truncnorm_fit")) stop("fit must be a HELP_truncnorm_fit object.")

  X_raw <- model.matrix(fit$formula, data = newdata)
  # Align columns defensively.
  missing_cols <- setdiff(names(fit$beta), colnames(X_raw))
  if (length(missing_cols) > 0L) stop("Missing columns in prediction matrix: ", paste(missing_cols, collapse = ", "))
  X_raw <- X_raw[, names(fit$beta), drop = FALSE]
  X <- sweep(X_raw, 2L, fit$x_center[names(fit$beta)], "-")
  X <- sweep(X, 2L, fit$x_scale[names(fit$beta)], "/")

  mu <- as.numeric(X %*% fit$beta)
  sigma <- fit$sigma
  z <- mu / sigma

  # E[Y | Y>0] = mu + sigma * phi(mu/sigma) / Phi(mu/sigma)
  log_mills <- dnorm(z, log = TRUE) - pnorm(z, log.p = TRUE)
  mills <- exp(log_mills)
  out <- mu + sigma * mills
  pmax(out, .Machine$double.eps)
}

HELP_aipw_details <- function(data, ps, ey1, ey0, conf_level = 0.95) {
  a <- as.numeric(data$treat)
  y <- as.numeric(data$dayslink)
  ps <- as.numeric(ps)
  ey1 <- as.numeric(ey1)
  ey0 <- as.numeric(ey0)

  if (!all(length(a) == c(length(y), length(ps), length(ey1), length(ey0)))) {
    stop("AIPW inputs must have the same length.")
  }

  # Numerical guard only; this does not intentionally trim practical positivity.
  eps <- 1e-10
  ps_safe <- pmin(pmax(ps, eps), 1 - eps)

  # EIF estimating function for E[Y(1)-Y(0)].
  psi_treated_minus_control <-
    (ey1 - ey0) +
    a / ps_safe * (y - ey1) -
    (1 - a) / (1 - ps_safe) * (y - ey0)

  # Paper estimand is control minus treated, so reverse the sign.
  contribution <- -psi_treated_minus_control
  estimate <- mean(contribution)
  se <- sd(contribution) / sqrt(length(contribution))
  zcrit <- qnorm(1 - (1 - conf_level) / 2)

  list(
    estimate = estimate,
    se = se,
    ci_low = estimate - zcrit * se,
    ci_high = estimate + zcrit * se,
    ps_min = min(ps),
    ps_max = max(ps),
    contribution = contribution
  )
}

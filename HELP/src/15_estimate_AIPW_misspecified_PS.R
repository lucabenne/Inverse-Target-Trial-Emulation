aipw_biased_prop = function(data, outcome_fit = NULL, return_details = FALSE){

  data_cov = data[,-c(1,3)]

  # Deliberately misspecified PS: omit three measured prognostic covariates.
  omit_ps = c("female", "daysdrink", "daysanysub")
  keep_ps = setdiff(names(data_cov), omit_ps)
  data_cov_biased = data_cov[, keep_ps, drop = FALSE]

  form_ps <- as.formula(paste("treat ~", paste(names(data_cov_biased), collapse = " + ")))
  ps <- as.numeric(predict(
    glm(form_ps, family = binomial, data = data),
    type = "response"
  ))

  # Correctly specified outcome nuisance: lower-truncated Normal regression.
  if (is.null(outcome_fit)) outcome_fit <- HELP_fit_truncnorm_outcome(data)

  ey1 <- HELP_predict_truncnorm_mean(
    outcome_fit,
    newdata = transform(data, treat = 1)
  )

  ey0 <- HELP_predict_truncnorm_mean(
    outcome_fit,
    newdata = transform(data, treat = 0)
  )

  out <- HELP_aipw_details(data, ps = ps, ey1 = ey1, ey0 = ey0)
  if (return_details) return(out)
  out$estimate
}

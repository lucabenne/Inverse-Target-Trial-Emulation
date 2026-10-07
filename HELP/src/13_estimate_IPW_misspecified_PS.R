ipw_biased = function(data,formula){

  data_cov = data[,-c(1,3)]

  # Deliberately misspecified PS: omit three measured prognostic covariates.
  # covariates used in the HELP simulation mechanisms.
  omit_ps = c("female", "daysdrink", "daysanysub")
  keep_ps = setdiff(names(data_cov), omit_ps)
  data_cov_biased = data_cov[, keep_ps, drop = FALSE]

  form_ps <- as.formula(paste("treat ~", paste(names(data_cov_biased), collapse = " + ")))
  ps <- as.numeric(predict(
    glm(form_ps, family = binomial, data = data),
    type = "response"
  ))

  eps <- 1e-10
  ps <- pmin(pmax(ps, eps), 1 - eps)

  ipw_value <- mean(
    data$treat * data$dayslink / ps -
      (1 - data$treat) * data$dayslink / (1 - ps)
  )

  -ipw_value
}

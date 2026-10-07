aipw = function(data, outcome_fit = NULL, return_details = FALSE){

  data_cov = data[,-c(1,3)]

  form_ps <- as.formula(paste("treat ~", paste(names(data_cov), collapse = " + ")))
  ps_fit <- glm(form_ps, family = binomial, data = data)
  ps <- as.numeric(predict(ps_fit, type = "response"))

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

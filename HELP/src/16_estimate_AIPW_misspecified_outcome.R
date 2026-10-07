aipw_biased_outcome = function(data, return_details = FALSE){

  data_cov = data[,-c(1,3)]
  data_cov_biased = data_cov[,-c(1,2,3), drop = FALSE]

  # Correct/full propensity-score model.
  form_ps <- as.formula(paste("treat ~", paste(names(data_cov), collapse = " + ")))
  ps <- as.numeric(predict(
    glm(form_ps, family = binomial, data = data),
    type = "response"
  ))

  # Deliberately misspecified outcome nuisance retained from the paper setup:
  # Gamma family with inverse link and omitted covariates.
  form_y <- as.formula(paste("dayslink ~ treat +", paste(names(data_cov_biased), collapse = " + ")))
  m <- glm(form_y, family = Gamma(link = "inverse"), data = data)

  ey1 <- as.numeric(predict(m, newdata = transform(data, treat = 1), type = "response"))
  ey0 <- as.numeric(predict(m, newdata = transform(data, treat = 0), type = "response"))

  out <- HELP_aipw_details(data, ps = ps, ey1 = ey1, ey0 = ey0)
  if (return_details) return(out)
  out$estimate
}

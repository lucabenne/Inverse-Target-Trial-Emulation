regression = function(data, outcome_fit = NULL){

  if (is.null(outcome_fit)) outcome_fit <- HELP_fit_truncnorm_outcome(data)

  ey1 <- HELP_predict_truncnorm_mean(
    outcome_fit,
    newdata = transform(data, treat = 1)
  )

  ey0 <- HELP_predict_truncnorm_mean(
    outcome_fit,
    newdata = transform(data, treat = 0)
  )

  # Paper estimand: control minus treated.
  mean(ey0 - ey1)
}

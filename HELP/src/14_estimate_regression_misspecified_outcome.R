regression_biased = function(data){
  
  data_cov = data[,-c(1,3)]
  data_cov_biased = data_cov[,-c(1,2,3)]
  
  
  form_y <- as.formula(paste("dayslink ~ treat +", paste(names(data_cov_biased), collapse = " + ")))
  m <- glm(form_y, family = Gamma(link = "inverse"), data = data)
  
  ey1 <- as.numeric(predict(
    m,
    newdata = transform(data, treat = 1),
    type = "response"
  ))
  
  ey0 <- as.numeric(predict(
    m,
    newdata = transform(data, treat = 0),
    type = "response"
  ))
  
  reg_value <- mean(ey1 - ey0)
  
  return(-reg_value)
  
}
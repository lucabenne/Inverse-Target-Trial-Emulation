ipw = function(data,formula){

library("ipw")
library("survey")
library("WeightIt")


data_cov = data[,-c(1,3)]
  

#compute ate through inverse probability weighting using package ipw


##########################################################################

form_ps <- as.formula(paste("treat ~", paste(names(data_cov), collapse = " + ")))

ps <- as.numeric(predict(
  glm(
    form_ps,
    family = binomial,
    data = data
  ),
  type = "response"
))

# --- IPW estimator ---
ipw <- function(a, y, ps) {
  mean(a * y / ps - (1 - a) * y / (1 - ps))
}

ipw_value <- ipw(
  a = data$treat,
  y = data$dayslink,
  ps = ps
)
ipw_value


return(-ipw_value)

}

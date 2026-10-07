overlapping_surv = function(data,index,sq=1){
  
  data_treat = data[which(data$treat==1),]
  data_nottreat = data[which(data$treat==0),]
  
  data_treat_cov = data_treat[,-index]
  data_nottreat_cov = data_nottreat[,-index]
  
  cov_treat = cov(data_treat_cov)
  cov_nottreat = cov(data_nottreat_cov)
  
  mu_treat = colMeans(data_treat_cov)
  mu_nottreat = colMeans(data_nottreat_cov)
  
  var_pop = solve((cov_treat+cov_nottreat)/2)
  
  overlapping = as.numeric(sqrt((mu_treat-mu_nottreat)%*%var_pop%*%(mu_treat-mu_nottreat)))
  overlapping2 = as.numeric((mu_treat-mu_nottreat)%*%var_pop%*%(mu_treat-mu_nottreat))
  
  if(sq==1){
  return(overlapping)
  }
  return(overlapping2)
}

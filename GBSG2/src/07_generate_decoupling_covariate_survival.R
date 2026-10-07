decoupling_covariate_surv1 = function(N,p,x,beta_weib=NULL,beta_weib_cens=NULL,sim_data){
  
  library(readr)
  library(MatchIt)
  library(readr)
  library(R2jags)
  library(truncnorm)
  library(sigmoid)
  

  
  data = sim_data
  
  
  #divide treated and untreated individuals
  
  
  ind_nottreat = which(data$treat==0)
  ind_treat = as.vector(which(data$treat==1))
  
  
  #order with respect to covariate x
  
  
  x_nottreat = as.data.frame(cbind(ind_nottreat,data[[x]][ind_nottreat]))
  x_nottreat = x_nottreat[order(x_nottreat$V2),]
  colnames(x_nottreat) = c("index","x")
  
  x_treat = as.data.frame(cbind(ind_treat,data[[x]][ind_treat]))
  x_treat = x_treat[order(-x_treat$V2),]
  colnames(x_treat) = c("index","x")
  
  
  #eliminate the first N untreated individuals (the ones with the lowest values of x) and
  #the last N treated individuals (the ones with the highest values of x)
  
  
  index_to_eliminate1 = x_nottreat$index[1:N]
  index_to_eliminate2 = x_treat$index[1:N]
  index_to_eliminate = c(index_to_eliminate1,index_to_eliminate2)
  
  
  data = data[-index_to_eliminate,]
  
  
  
  return(data)
}

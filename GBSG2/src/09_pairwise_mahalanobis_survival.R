pairwise_mahalanobis = function(datt,index){
  
  data1_cov = datt[,-index]
  
  ncol = ncol(data1_cov)
  
  ind_nottreat = which(datt$treat==0)
  ind_treat = as.vector(which(datt$treat==1))
  
  treat = data1_cov[ind_treat,]
  nottreat = data1_cov[ind_nottreat,]
  
  length_nottreat = length(ind_nottreat)
  length_treat = length(ind_treat)
  
  d = matrix(0,ncol,length_nottreat)
  Ma = matrix(0,length_nottreat,length_treat)
  
  for (j in 1:length_treat) {
    
    
    for (i in 1:length_nottreat){
      
      d[,i] = as.matrix(as.numeric(treat[j,])-as.numeric(nottreat[i,]))[,1]
      
      e = solve(cov(data1_cov),d[,i])
      
      Ma[i,j] = sum(d[,i]*e)
    }
  }
  
  return(Ma)
}
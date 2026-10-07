library(MatchIt)
library(marginaleffects)

tte_mahalanobis = function(data, formula){

  
#using MatchIt package compute pairs of individuals (starting from the first treated one) from different populations
#with the minimum mahalanobis distance  
   
m.out = matchit(formula=formula,
                   data=data,
                method = "full",
                   distance='mahalanobis')

m.data <- match.data(m.out, data = data)

head(m.data)

fit <- lm(dayslink ~ treat + female + daysdrink + daysanysub + pcs + mcs + 
                            cesd + sexrisk + drugrisk, data = m.data, weights = weights)

a=avg_comparisons(fit,
                variables = "treat",
                vcov = ~subclass,
                newdata = m.data,
                wts = "weights")


return(-as.numeric(a[3]))

#nrow = nrow(data)
#a=1:nrow


#construct a dataframe x with the indexes of coupled pairs


#x = m.out$match.matrix
#x = as.numeric(x[,1])
#y=m.out$treat
#y = cbind(y,1:nrow)
#y=as.data.frame(y)
#y = subset(y, y==1)
#z = cbind(y,x)
#x = z[,2:3]
#x = as.data.frame(x)


#select uncoupled individuals in case there are more untreated individuals with respect to the treated ones


#if (sum(data$treat) >= nrow-sum(data$treat)){
#not_matched = x$V2[which(is.na(x$x)==TRUE)]
#}


#select uncoupled individuals in case there are more treated individuals with respect to the untreated ones

  
#if (sum(data$treat) < nrow-sum(data$treat)){
#  b=c(x$V2,x$x)
#  not_matched = a[-b]
#}


#eliminate uncoupled individuals from the data and compute att


#tte = data[-not_matched,]
#tte_treat = tte[which(tte$treat==1),]
#tte_nottreat = tte[which(tte$treat==0),]
#att = mean(tte_nottreat$dayslink - tte_treat$dayslink)

#data_att = list(tte, att)

#return(data_att)

}


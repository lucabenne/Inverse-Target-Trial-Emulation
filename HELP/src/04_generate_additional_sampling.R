if (!exists("HELP_generate_outcome")) source("03a_generate_HELP_synthetic_trial.R")

newdatamahal_normal2 = function(data, bern, n, q, BG){

  library(readr)
  library(MatchIt)
  library(truncnorm)

  data_cov = data[,-c(1,3)]
  data_mahal = data[,-c(3)]
  data_cov_out = data[,-c(1)]

  help_clean = data
  help_clean_cov = data_cov

  # Match treated and untreated individuals in the source HELP RCT.
  formula <- as.formula("treat ~ female + daysdrink + daysanysub + pcs + mcs + cesd + sexrisk + drugrisk")

  m.out = matchit(
    formula = formula,
    data = data_mahal,
    distance = 'mahalanobis'
  )

  nrow_data = nrow(data)

  x = m.out$match.matrix
  x = as.numeric(x[,1])
  y = m.out$treat
  y = cbind(y,1:nrow_data)
  y = as.data.frame(y)
  y = subset(y, y==1)
  z = cbind(y,x)
  x = z[,2:3]
  x = as.data.frame(x)

  ncol_cov = ncol(help_clean_cov)

  ind_nottreat = which(data$treat==0)
  ind_treat = as.vector(which(data$treat==1))

  treat_cov = help_clean_cov[ind_treat,]
  nottreat_cov = help_clean_cov[ind_nottreat,]

  length_nottreat = length(ind_nottreat)
  length_treat = length(ind_treat)

  # Pairwise Mahalanobis distance between treated and untreated source subjects.
  d = matrix(0,ncol_cov,length_nottreat)
  Ma = matrix(0,length_nottreat,length_treat)

  for (jj in 1:length_treat) {
    for (ii in 1:length_nottreat){
      d[,ii] = as.matrix(as.numeric(treat_cov[jj,])-as.numeric(nottreat_cov[ii,]))[,1]
      e = solve(cov(help_clean_cov),d[,ii])
      Ma[ii,jj] = sum(d[,ii]*e)
    }
  }

  list_data = list()

  for (j in 1:length(n)) {

    # The Additional Sampling mechanism perturbs source-RCT covariates locally.
    # Outcome-model parameters are fixed at posterior means so the treatment-effect
    # target remains constant across Monte Carlo replicates.

    Ma_ordered = sort(Ma)
    quantile_value = Ma_ordered[round(length(Ma_ordered)*q[j])]

    selected_pair = which(Ma == quantile_value, arr.ind = TRUE)[1, ]
    row = selected_pair[1]
    colmn = selected_pair[2]

    max_treated = data_cov_out[ind_treat[colmn],]
    max_nottreated = data_cov_out[ind_nottreat[row],]

    data = help_clean

    for (i in 1:n[j]) {

      data_cov_current = data[,-c(1)]
      mu = as.data.frame(t(colMeans(data_cov_current)))
      next_treat = mu
      next_nottreat = mu

      next_treat$female = rbinom(1,1,(max_treated$female+0.5)/2)
      next_treat$daysdrink = round(rtruncnorm(1,0,Inf,max_treated$daysdrink,
                                               abs(max_treated$daysdrink-mu$daysdrink)/4))
      next_treat$daysanysub = round(rtruncnorm(1,0,Inf,max_treated$daysanysub,
                                                abs(max_treated$daysanysub-mu$daysanysub)/4))
      next_treat$cesd = round(rtruncnorm(1,0,Inf,max_treated$cesd,
                                          abs(max_treated$cesd-mu$cesd)/4))
      next_treat$sexrisk = round(rtruncnorm(1,0,Inf,max_treated$sexrisk,
                                             abs(max_treated$sexrisk-mu$sexrisk)/4))
      next_treat$drugrisk = round(rtruncnorm(1,0,Inf,max_treated$drugrisk,
                                              abs(max_treated$drugrisk-mu$drugrisk)/4))
      next_treat$pcs = rtruncnorm(1,0,Inf,max_treated$pcs,
                                   abs(max_treated$pcs-mu$pcs)/4)
      next_treat$mcs = rtruncnorm(1,0,Inf,max_treated$mcs,
                                   abs(max_treated$mcs-mu$mcs)/4)

      cov_treat = next_treat[,c("female","daysdrink","daysanysub","pcs","mcs",
                                 "cesd","sexrisk","drugrisk"),drop=FALSE]
      next_treat$dayslink = HELP_generate_outcome(
        covariates = cov_treat,
        treat = bern[i],
        BG = BG,
        posterior_mode = "mean",
        posterior_draw = NULL,
        round_outcome = TRUE
      )

      next_nottreat$female = rbinom(1,1,(max_nottreated$female+0.5)/2)
      next_nottreat$daysdrink = round(rtruncnorm(1,0,Inf,max_nottreated$daysdrink,
                                                  abs(max_nottreated$daysdrink-mu$daysdrink)/4))
      next_nottreat$daysanysub = round(rtruncnorm(1,0,Inf,max_nottreated$daysanysub,
                                                   abs(max_nottreated$daysanysub-mu$daysanysub)/4))
      next_nottreat$cesd = round(rtruncnorm(1,0,Inf,max_nottreated$cesd,
                                             abs(max_nottreated$cesd-mu$cesd)/4))
      next_nottreat$sexrisk = round(rtruncnorm(1,0,Inf,max_nottreated$sexrisk,
                                                abs(max_nottreated$sexrisk-mu$sexrisk)/4))
      next_nottreat$drugrisk = round(rtruncnorm(1,0,Inf,max_nottreated$drugrisk,
                                                 abs(max_nottreated$drugrisk-mu$drugrisk)/4))
      next_nottreat$pcs = rtruncnorm(1,0,Inf,max_nottreated$pcs,
                                      abs(max_nottreated$pcs-mu$pcs)/4)
      next_nottreat$mcs = rtruncnorm(1,0,Inf,max_nottreated$mcs,
                                      abs(max_nottreated$mcs-mu$mcs)/4)

      cov_nottreat = next_nottreat[,c("female","daysdrink","daysanysub","pcs","mcs",
                                       "cesd","sexrisk","drugrisk"),drop=FALSE]
      next_nottreat$dayslink = HELP_generate_outcome(
        covariates = cov_nottreat,
        treat = bern[i],
        BG = BG,
        posterior_mode = "mean",
        posterior_draw = NULL,
        round_outcome = TRUE
      )

      if(bern[i] == 1){
        next_treat_comp = data.frame(
          treat = 1,
          female = next_treat$female,
          dayslink = next_treat$dayslink,
          daysdrink = next_treat$daysdrink,
          daysanysub = next_treat$daysanysub,
          pcs = next_treat$pcs,
          mcs = next_treat$mcs,
          cesd = next_treat$cesd,
          sexrisk = next_treat$sexrisk,
          drugrisk = next_treat$drugrisk
        )
        data = rbind(data,next_treat_comp)
      }

      if(bern[i] == 0){
        next_nottreat_comp = data.frame(
          treat = 0,
          female = next_nottreat$female,
          dayslink = next_nottreat$dayslink,
          daysdrink = next_nottreat$daysdrink,
          daysanysub = next_nottreat$daysanysub,
          pcs = next_nottreat$pcs,
          mcs = next_nottreat$mcs,
          cesd = next_nottreat$cesd,
          sexrisk = next_nottreat$sexrisk,
          drugrisk = next_nottreat$drugrisk
        )
        data = rbind(data,next_nottreat_comp)
      }
    }

    rownames(data) = NULL
    list_data = append(list_data,list(data))
  }

  return(list_data)
}

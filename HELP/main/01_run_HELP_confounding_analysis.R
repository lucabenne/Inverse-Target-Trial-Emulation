library(readr)
library(MatchIt)
library(readr)
library(R2jags)
library(truncnorm)
library(sigmoid)
#library(rMR)
library(ggplot2)
library(xtable)

source("../src/09_estimate_mahalanobis_matching.R")
source("../src/06_generate_perfect_knowledge.R")
source("../src/08_measure_mahalanobis_imbalance.R")
source("../src/05_generate_decoupling_mahalanobis.R")
source("../src/07_generate_decoupling_covariate.R")
source("../src/04_generate_additional_sampling.R")
source("../src/10_estimate_IPW.R")
source("../src/11a_fit_truncated_normal_outcome.R")
source("../src/12_estimate_AIPW.R")
source("../src/11_estimate_outcome_regression.R")
source("../src/02_fit_HELP_bayesian_DGP.R")
source("../src/03a_generate_HELP_synthetic_trial.R")
source("../src/13_estimate_IPW_misspecified_PS.R")
source("../src/14_estimate_regression_misspecified_outcome.R")
source("../src/15_estimate_AIPW_misspecified_PS.R")
source("../src/16_estimate_AIPW_misspecified_outcome.R")





formula = as.formula("treat ~ female + daysdrink + daysanysub + pcs + mcs + cesd + sexrisk + drugrisk")

# Monte Carlo replicates used in the paper
n_rep = 1000L

help_clean = read_csv("../HELP_analysis_data.csv", show_col_types = FALSE)
analysis_vars = c("treat","female","dayslink","daysdrink","daysanysub","pcs","mcs","cesd","sexrisk","drugrisk")
data = as.data.frame(help_clean[,analysis_vars])


bern = rbinom(n_rep,1,0.66)

iseed = 101010L
set.seed(iseed)
if (!exists("BG", inherits = FALSE) || is.null(BG$HELP_source_n) || BG$HELP_source_n != nrow(data)) {
  BG = jags_data(data, seed = iseed)
} else {
  cat("\nUsing pre-fitted HELP Bayesian DGP from this session.\n")
}

# Reset so downstream simulation starts from the intended seed independently of JAGS RNG use.
set.seed(iseed)


overlapping0 = overlapping(data,3)
# Source-RCT benchmark used throughout the paper analyses.
# theta91 is a latent-location coefficient under the lower-truncated outcome DGP
# and is therefore not itself the observed-scale mean treatment effect.
ate_before0 = mean(data[data$treat==0,]$dayslink) - mean(data[data$treat==1,]$dayslink)
att_tte0 = tte_mahalanobis(data, formula)
ate_ipw0 = ipw(data, formula)
outcome_fit0 = HELP_fit_truncnorm_outcome(data)
ate_reg0 = regression(data, outcome_fit = outcome_fit0)
aipw0_details = aipw(data, outcome_fit = outcome_fit0, return_details = TRUE)
ate_aipw0 = aipw0_details$estimate

ate0 = c(ate_before0,att_tte0,ate_ipw0, ate_reg0, ate_aipw0)



iseed = 101010
set.seed(iseed)

overlapping2=numeric()
ate_before2=numeric()
att_tte2=numeric()
ate_ipw2=numeric()
ate_reg2=numeric()
ate_aipw2=numeric()
aipw_se2=numeric()
aipw_ci_low2=numeric()
aipw_ci_high2=numeric()
outcome_fit2=vector("list", n_rep)
sum_abs_error2 = numeric()

# Final Additional Sampling calibration used for the revised manuscript.
q = rep(0.90, n_rep)
# Add 60 individuals to the source HELP RCT in every Monte Carlo replicate.
n_newdatamahal = rep(60L, n_rep)


list_data2 = newdatamahal_normal2(data,rbinom(n_rep,1,0.6),n_newdatamahal,q, BG)
# Additional Sampling augments the actual source RCT, so its benchmark stays fixed.
target2 = rep(ate_before0, n_rep)

for (i in seq_len(n_rep)) {
data2 = list_data2[[i]]  

overlapping2[i] = overlapping(data2,3)
ate_before2[i] = mean(data2[which(data2$treat==0),]$dayslink) - mean(data2[which(data2$treat==1),]$dayslink)
att_tte2[i] = tte_mahalanobis(data2,formula)
ate_ipw2[i] = ipw(data2, formula)
outcome_fit2[[i]] = HELP_fit_truncnorm_outcome(data2)
ate_reg2[i] = regression(data2, outcome_fit = outcome_fit2[[i]])
aipw_details = aipw(data2, outcome_fit = outcome_fit2[[i]], return_details = TRUE)
ate_aipw2[i] = aipw_details$estimate
aipw_se2[i] = aipw_details$se
aipw_ci_low2[i] = aipw_details$ci_low
aipw_ci_high2[i] = aipw_details$ci_high
print(i)
}

ate2 = c(mean(ate_before2),mean(att_tte2),mean(ate_ipw2), mean(ate_reg2), mean(ate_aipw2))
var_ate2 = c(var(ate_before2),var(att_tte2),var(ate_ipw2),var(ate_reg2),var(ate_aipw2))

sum_abs_error2[1] = mean(abs(ate_before2-target2))
sum_abs_error2[2] = mean(abs(att_tte2-target2))
sum_abs_error2[3] = mean(abs(ate_ipw2-target2))
sum_abs_error2[4] = mean(abs(ate_reg2-target2))
sum_abs_error2[5] = mean(abs(ate_aipw2-target2))


iseed = 101010
set.seed(iseed)

overlapping3=numeric()
ate_before3=numeric()
att_tte3=numeric()
ate_ipw3=numeric()
ate_reg3=numeric()
ate_aipw3=numeric()
aipw_se3=numeric()
aipw_ci_low3=numeric()
aipw_ci_high3=numeric()
outcome_fit3=vector("list", n_rep)
sum_abs_error3 = numeric()

resampling=0
n_res = 1000

# Moderate Decoupling Mahalanobis scenario.
N = rep(150L,n_rep)

list_data3 = decoupling_mahalanobis2(data,N,0.54,'daysanysub',resampling,n_res, BG)
# Each decoupling replicate is paired with the randomized synthetic RCT from
# which it was created, exactly as in the GBSG2 decoupling analysis.
target3 = vapply(list_data3, function(d) attr(d, "rct_target"), numeric(1))
if (any(!is.finite(target3))) stop("Missing Decoupling Mahalanobis replicate-specific RCT targets.")



for (i in seq_len(n_rep)) {

data3 = list_data3[[i]]

overlapping3[i] = overlapping(data3,3)
ate_before3[i] = mean(data3[which(data3$treat==0),]$dayslink) - mean(data3[which(data3$treat==1),]$dayslink)
att_tte3[i] = tte_mahalanobis(data3,formula)
ate_ipw3[i] = ipw(data3, formula)
outcome_fit3[[i]] = HELP_fit_truncnorm_outcome(data3)
ate_reg3[i] = regression(data3, outcome_fit = outcome_fit3[[i]])
aipw_details = aipw(data3, outcome_fit = outcome_fit3[[i]], return_details = TRUE)
ate_aipw3[i] = aipw_details$estimate
aipw_se3[i] = aipw_details$se
aipw_ci_low3[i] = aipw_details$ci_low
aipw_ci_high3[i] = aipw_details$ci_high

print(i)
}

ate3 = c(mean(ate_before3),mean(att_tte3),mean(ate_ipw3), mean(ate_reg3), mean(ate_aipw3))
var_ate3 = c(var(ate_before3),var(att_tte3),var(ate_ipw3),var(ate_reg3),var(ate_aipw3))


sum_abs_error3[1] = mean(abs(ate_before3-target3))
sum_abs_error3[2] = mean(abs(att_tte3-target3))
sum_abs_error3[3] = mean(abs(ate_ipw3-target3))
sum_abs_error3[4] = mean(abs(ate_reg3-target3))
sum_abs_error3[5] = mean(abs(ate_aipw3-target3))



iseed = 101010
set.seed(iseed)

overlapping4=numeric()
ate_before4=numeric()
att_tte4=numeric()
ate_ipw4=numeric()
ate_reg4=numeric()
ate_aipw4=numeric()
aipw_se4=numeric()
aipw_ci_low4=numeric()
aipw_ci_high4=numeric()
outcome_fit4=vector("list", n_rep)
sum_abs_error4 = numeric()

# Moderate Perfect Knowledge scenario: fixed sample size and reduced
# treatment-allocation strength to preserve meaningful overlap.
n = rep(1000L,n_rep)
perfect_gamma = 1.25

list_data4 = observational_help2(data=data, n=n, BG=BG, gamma=perfect_gamma)
# Perfect Knowledge gets an auxiliary randomized synthetic-RCT benchmark on
# the same generated covariate population before treatment is made dependent on X.
target4 = vapply(list_data4, function(d) attr(d, "rct_target"), numeric(1))
if (any(!is.finite(target4))) stop("Missing Perfect Knowledge replicate-specific RCT targets.")

for (i in seq_len(n_rep)) {

data4 = list_data4[[i]]
  
overlapping4[i] = overlapping(data4,3)
ate_before4[i] = mean(data4[which(data4$treat==0),]$dayslink) - mean(data4[which(data4$treat==1),]$dayslink)
att_tte4[i] = tte_mahalanobis(data4,formula)
ate_ipw4[i] = ipw(data4, formula)
outcome_fit4[[i]] = HELP_fit_truncnorm_outcome(data4)
ate_reg4[i] = regression(data4, outcome_fit = outcome_fit4[[i]])
aipw_details = aipw(data4, outcome_fit = outcome_fit4[[i]], return_details = TRUE)
ate_aipw4[i] = aipw_details$estimate
aipw_se4[i] = aipw_details$se
aipw_ci_low4[i] = aipw_details$ci_low
aipw_ci_high4[i] = aipw_details$ci_high

print(i)
}


ate4 = c(mean(ate_before4),mean(att_tte4),mean(ate_ipw4), mean(ate_reg4), mean(ate_aipw4))
var_ate4 = c(var(ate_before4),var(att_tte4),var(ate_ipw4),var(ate_reg4),var(ate_aipw4))


sum_abs_error4[1] = mean(abs(ate_before4-target4))
sum_abs_error4[2] = mean(abs(att_tte4-target4))
sum_abs_error4[3] = mean(abs(ate_ipw4-target4))
sum_abs_error4[4] = mean(abs(ate_reg4-target4))
sum_abs_error4[5] = mean(abs(ate_aipw4-target4))



iseed = 101010
set.seed(iseed)

overlapping5=numeric()
ate_before5=numeric()
att_tte5=numeric()
ate_ipw5=numeric()
ate_reg5=numeric()
ate_aipw5=numeric()
aipw_se5=numeric()
aipw_ci_low5=numeric()
aipw_ci_high5=numeric()
outcome_fit5=vector("list", n_rep)
sum_abs_error5 = numeric()

# Moderate Decoupling Covariate scenario.
N = rep(90L,n_rep)

list_data5 = decoupling_covariate(data,N,0.5,'daysanysub', BG)
# Paired randomized synthetic-RCT benchmark before covariate decoupling.
target5 = vapply(list_data5, function(d) attr(d, "rct_target"), numeric(1))
if (any(!is.finite(target5))) stop("Missing Decoupling Covariate replicate-specific RCT targets.")


for (i in seq_len(n_rep)) {

data5 = list_data5[[i]]
  
overlapping5[i] = overlapping(data5,3)
ate_before5[i] = mean(data5[which(data5$treat==0),]$dayslink) - mean(data5[which(data5$treat==1),]$dayslink)
att_tte5[i] = tte_mahalanobis(data5,formula)
ate_ipw5[i] = ipw(data5, formula)
outcome_fit5[[i]] = HELP_fit_truncnorm_outcome(data5)
ate_reg5[i] = regression(data5, outcome_fit = outcome_fit5[[i]])
aipw_details = aipw(data5, outcome_fit = outcome_fit5[[i]], return_details = TRUE)
ate_aipw5[i] = aipw_details$estimate
aipw_se5[i] = aipw_details$se
aipw_ci_low5[i] = aipw_details$ci_low
aipw_ci_high5[i] = aipw_details$ci_high

print(i)
}


ate5 = c(mean(ate_before5),mean(att_tte5),mean(ate_ipw5), mean(ate_reg5), mean(ate_aipw5))
var_ate5 = c(var(ate_before5),var(att_tte5),var(ate_ipw5),var(ate_reg5),var(ate_aipw5))


sum_abs_error5[1] = mean(abs(ate_before5-target5))
sum_abs_error5[2] = mean(abs(att_tte5-target5))
sum_abs_error5[3] = mean(abs(ate_ipw5-target5))
sum_abs_error5[4] = mean(abs(ate_reg5-target5))
sum_abs_error5[5] = mean(abs(ate_aipw5-target5))


ate = data.frame(cbind(ate0,ate2,ate3,ate4,ate5))
rownames(ate) = c("ATE_before","ATE_Mahalanobis","ATE_ipw","ATE_regression", "ATE_aipw")
colnames(ate) = c("Original RCT", "Rand new dataset", "Dec Mahalanobis", "Perfect Knowledge", "Dec covariate")

sum_abs_error = data.frame(cbind(sum_abs_error2,sum_abs_error3,sum_abs_error4,sum_abs_error5))
colnames(sum_abs_error) = c("Randomised new dataset", "Decoupling Mahalanobis", "Perfect Knowledge", "Decoupling covariate")
rownames(sum_abs_error) = c("ATE_before","ATE_Mahalanobis","ATE_ipw","ATE_regression", "ATE_aipw")


######################################################################

#double robustness property AIPW

######################################################################


ate_before0_biased_outcome = ate_before0


iseed = 101010
set.seed(iseed)


ate_before2_biased_outcome=numeric()
ate_reg2_biased_outcome=numeric()
ate_aipw2_biased_outcome=numeric()
aipw_se2_biased_outcome=numeric()
aipw_ci_low2_biased_outcome=numeric()
aipw_ci_high2_biased_outcome=numeric()
sum_abs_error2_biased_outcome = numeric()



for (i in seq_len(n_rep)) {
  data2 = list_data2[[i]]  
  
  ate_before2_biased_outcome[i] = mean(data2[which(data2$treat==0),]$dayslink) - mean(data2[which(data2$treat==1),]$dayslink)
  ate_reg2_biased_outcome[i] = regression_biased(data2)
  aipw_details = aipw_biased_outcome(data2, return_details = TRUE)
  ate_aipw2_biased_outcome[i] = aipw_details$estimate
  aipw_se2_biased_outcome[i] = aipw_details$se
  aipw_ci_low2_biased_outcome[i] = aipw_details$ci_low
  aipw_ci_high2_biased_outcome[i] = aipw_details$ci_high
  print(i)
}

ate2_biased_outcome = c(mean(ate_before2_biased_outcome),mean(ate_reg2_biased_outcome), mean(ate_aipw2_biased_outcome))
var_ate2_biased_outcome = c(var(ate_before2_biased_outcome),var(ate_reg2_biased_outcome),var(ate_aipw2_biased_outcome))

sum_abs_error2_biased_outcome[1] = mean(abs(ate_before2_biased_outcome-target2))
sum_abs_error2_biased_outcome[2] = mean(abs(ate_reg2_biased_outcome-target2))
sum_abs_error2_biased_outcome[3] = mean(abs(ate_aipw2_biased_outcome-target2))


iseed = 101010
set.seed(iseed)


ate_before3_biased_outcome=numeric()
ate_reg3_biased_outcome=numeric()
ate_aipw3_biased_outcome=numeric()
aipw_se3_biased_outcome=numeric()
aipw_ci_low3_biased_outcome=numeric()
aipw_ci_high3_biased_outcome=numeric()
sum_abs_error3_biased_outcome = numeric()



for (i in seq_len(n_rep)) {
  
  data3 = list_data3[[i]]
  
  ate_before3_biased_outcome[i] = mean(data3[which(data3$treat==0),]$dayslink) - mean(data3[which(data3$treat==1),]$dayslink)
  ate_reg3_biased_outcome[i] = regression_biased(data3)
  aipw_details = aipw_biased_outcome(data3, return_details = TRUE)
  ate_aipw3_biased_outcome[i] = aipw_details$estimate
  aipw_se3_biased_outcome[i] = aipw_details$se
  aipw_ci_low3_biased_outcome[i] = aipw_details$ci_low
  aipw_ci_high3_biased_outcome[i] = aipw_details$ci_high
  
  print(i)
}

ate3_biased_outcome = c(mean(ate_before3_biased_outcome), mean(ate_reg3_biased_outcome), mean(ate_aipw3_biased_outcome))
var_ate3_biased_outcome = c(var(ate_before3_biased_outcome),var(ate_reg3_biased_outcome),var(ate_aipw3_biased_outcome))


sum_abs_error3_biased_outcome[1] = mean(abs(ate_before3_biased_outcome-target3))
sum_abs_error3_biased_outcome[2] = mean(abs(ate_reg3_biased_outcome-target3))
sum_abs_error3_biased_outcome[3] = mean(abs(ate_aipw3_biased_outcome-target3))



iseed = 101010
set.seed(iseed)

ate_before4_biased_outcome=numeric()
ate_reg4_biased_outcome=numeric()
ate_aipw4_biased_outcome=numeric()
aipw_se4_biased_outcome=numeric()
aipw_ci_low4_biased_outcome=numeric()
aipw_ci_high4_biased_outcome=numeric()
sum_abs_error4_biased_outcome = numeric()



for (i in seq_len(n_rep)) {
  
  data4 = list_data4[[i]]
  
  ate_before4_biased_outcome[i] = mean(data4[which(data4$treat==0),]$dayslink) - mean(data4[which(data4$treat==1),]$dayslink)
  ate_reg4_biased_outcome[i] = regression_biased(data4)
  aipw_details = aipw_biased_outcome(data4, return_details = TRUE)
  ate_aipw4_biased_outcome[i] = aipw_details$estimate
  aipw_se4_biased_outcome[i] = aipw_details$se
  aipw_ci_low4_biased_outcome[i] = aipw_details$ci_low
  aipw_ci_high4_biased_outcome[i] = aipw_details$ci_high
  
  print(i)
}


ate4_biased_outcome = c(mean(ate_before4_biased_outcome),mean(ate_reg4_biased_outcome), mean(ate_aipw4_biased_outcome))
var_ate4_biased_outcome = c(var(ate_before4_biased_outcome),var(ate_reg4_biased_outcome),var(ate_aipw4_biased_outcome))


sum_abs_error4_biased_outcome[1] = mean(abs(ate_before4_biased_outcome-target4))
sum_abs_error4_biased_outcome[2] = mean(abs(ate_reg4_biased_outcome-target4))
sum_abs_error4_biased_outcome[3] = mean(abs(ate_aipw4_biased_outcome-target4))



iseed = 101010
set.seed(iseed)

ate_before5_biased_outcome=numeric()
ate_reg5_biased_outcome=numeric()
ate_aipw5_biased_outcome=numeric()
aipw_se5_biased_outcome=numeric()
aipw_ci_low5_biased_outcome=numeric()
aipw_ci_high5_biased_outcome=numeric()
sum_abs_error5_biased_outcome = numeric()



for (i in seq_len(n_rep)) {
  
  data5 = list_data5[[i]]
  
  ate_before5_biased_outcome[i] = mean(data5[which(data5$treat==0),]$dayslink) - mean(data5[which(data5$treat==1),]$dayslink)
  ate_reg5_biased_outcome[i] = regression_biased(data5)
  aipw_details = aipw_biased_outcome(data5, return_details = TRUE)
  ate_aipw5_biased_outcome[i] = aipw_details$estimate
  aipw_se5_biased_outcome[i] = aipw_details$se
  aipw_ci_low5_biased_outcome[i] = aipw_details$ci_low
  aipw_ci_high5_biased_outcome[i] = aipw_details$ci_high
  
  print(i)
}


ate5_biased_outcome = c(mean(ate_before5_biased_outcome),mean(ate_reg5_biased_outcome), mean(ate_aipw5_biased_outcome))
var_ate5_biased_outcome = c(var(ate_before5_biased_outcome),var(ate_reg5_biased_outcome),var(ate_aipw5_biased_outcome))


sum_abs_error5_biased_outcome[1] = mean(abs(ate_before5_biased_outcome-target5))
sum_abs_error5_biased_outcome[2] = mean(abs(ate_reg5_biased_outcome-target5))
sum_abs_error5_biased_outcome[3] = mean(abs(ate_aipw5_biased_outcome-target5))



# Monte Carlo simulation intervals (diagnostic only)

as.numeric(quantile(ate_reg2_biased_outcome,c(0.025,0.975)))
as.numeric(quantile(ate_aipw2_biased_outcome,c(0.025,0.975)))

as.numeric(quantile(ate_reg3_biased_outcome,c(0.025,0.975)))
as.numeric(quantile(ate_aipw3_biased_outcome,c(0.025,0.975)))

as.numeric(quantile(ate_reg4_biased_outcome,c(0.025,0.975)))
as.numeric(quantile(ate_aipw4_biased_outcome,c(0.025,0.975)))

as.numeric(quantile(ate_reg5_biased_outcome,c(0.025,0.975)))
as.numeric(quantile(ate_aipw5_biased_outcome,c(0.025,0.975)))



sum_abs_error5_biased_outcome[1] = mean(abs(ate_before5_biased_outcome-target5))
sum_abs_error5_biased_outcome[2] = mean(abs(ate_reg5_biased_outcome-target5))
sum_abs_error5_biased_outcome[3] = mean(abs(ate_aipw5_biased_outcome-target5))



#################################################################################

#biased propensity score

#################################################################################



ate_before0_biased_prop = ate_before0


iseed = 101010
set.seed(iseed)


ate_before2_biased_prop=numeric()
ate_ipw2_biased_prop=numeric()
ate_aipw2_biased_prop=numeric()
aipw_se2_biased_prop=numeric()
aipw_ci_low2_biased_prop=numeric()
aipw_ci_high2_biased_prop=numeric()
sum_abs_error2_biased_prop = numeric()



for (i in seq_len(n_rep)) {
  data2 = list_data2[[i]]  
  
  ate_before2_biased_prop[i] = mean(data2[which(data2$treat==0),]$dayslink) - mean(data2[which(data2$treat==1),]$dayslink)
  ate_ipw2_biased_prop[i] = ipw_biased(data2,formula)
  aipw_details = aipw_biased_prop(data2, outcome_fit = outcome_fit2[[i]], return_details = TRUE)
  ate_aipw2_biased_prop[i] = aipw_details$estimate
  aipw_se2_biased_prop[i] = aipw_details$se
  aipw_ci_low2_biased_prop[i] = aipw_details$ci_low
  aipw_ci_high2_biased_prop[i] = aipw_details$ci_high
  print(i)
}

ate2_biased_prop = c(mean(ate_before2_biased_prop),mean(ate_ipw2_biased_prop), mean(ate_aipw2_biased_prop))
var_ate2_biased_prop = c(var(ate_before2_biased_prop),var(ate_ipw2_biased_prop),var(ate_aipw2_biased_prop))

sum_abs_error2_biased_prop[1] = mean(abs(ate_before2_biased_prop-target2))
sum_abs_error2_biased_prop[2] = mean(abs(ate_ipw2_biased_prop-target2))
sum_abs_error2_biased_prop[3] = mean(abs(ate_aipw2_biased_prop-target2))


iseed = 101010
set.seed(iseed)


ate_before3_biased_prop=numeric()
ate_ipw3_biased_prop=numeric()
ate_aipw3_biased_prop=numeric()
aipw_se3_biased_prop=numeric()
aipw_ci_low3_biased_prop=numeric()
aipw_ci_high3_biased_prop=numeric()
sum_abs_error3_biased_prop = numeric()



for (i in seq_len(n_rep)) {
  
  data3 = list_data3[[i]]
  
  ate_before3_biased_prop[i] = mean(data3[which(data3$treat==0),]$dayslink) - mean(data3[which(data3$treat==1),]$dayslink)
  ate_ipw3_biased_prop[i] = ipw_biased(data3,formula)
  aipw_details = aipw_biased_prop(data3, outcome_fit = outcome_fit3[[i]], return_details = TRUE)
  ate_aipw3_biased_prop[i] = aipw_details$estimate
  aipw_se3_biased_prop[i] = aipw_details$se
  aipw_ci_low3_biased_prop[i] = aipw_details$ci_low
  aipw_ci_high3_biased_prop[i] = aipw_details$ci_high
  
  print(i)
}

ate3_biased_prop = c(mean(ate_before3_biased_prop), mean(ate_ipw3_biased_prop), mean(ate_aipw3_biased_prop))
var_ate3_biased_prop = c(var(ate_before3_biased_prop),var(ate_ipw3_biased_prop),var(ate_aipw3_biased_prop))


sum_abs_error3_biased_prop[1] = mean(abs(ate_before3_biased_prop-target3))
sum_abs_error3_biased_prop[2] = mean(abs(ate_ipw3_biased_prop-target3))
sum_abs_error3_biased_prop[3] = mean(abs(ate_aipw3_biased_prop-target3))



iseed = 101010
set.seed(iseed)

ate_before4_biased_prop=numeric()
ate_ipw4_biased_prop=numeric()
ate_aipw4_biased_prop=numeric()
aipw_se4_biased_prop=numeric()
aipw_ci_low4_biased_prop=numeric()
aipw_ci_high4_biased_prop=numeric()
sum_abs_error4_biased_prop = numeric()



for (i in seq_len(n_rep)) {
  
  data4 = list_data4[[i]]
  
  ate_before4_biased_prop[i] = mean(data4[which(data4$treat==0),]$dayslink) - mean(data4[which(data4$treat==1),]$dayslink)
  ate_ipw4_biased_prop[i] = ipw_biased(data4,formula)
  aipw_details = aipw_biased_prop(data4, outcome_fit = outcome_fit4[[i]], return_details = TRUE)
  ate_aipw4_biased_prop[i] = aipw_details$estimate
  aipw_se4_biased_prop[i] = aipw_details$se
  aipw_ci_low4_biased_prop[i] = aipw_details$ci_low
  aipw_ci_high4_biased_prop[i] = aipw_details$ci_high
  
  print(i)
}


ate4_biased_prop = c(mean(ate_before4_biased_prop),mean(ate_ipw4_biased_prop), mean(ate_aipw4_biased_prop))
var_ate4_biased_prop = c(var(ate_before4_biased_prop),var(ate_ipw4_biased_prop),var(ate_aipw4_biased_prop))


sum_abs_error4_biased_prop[1] = mean(abs(ate_before4_biased_prop-target4))
sum_abs_error4_biased_prop[2] = mean(abs(ate_ipw4_biased_prop-target4))
sum_abs_error4_biased_prop[3] = mean(abs(ate_aipw4_biased_prop-target4))



iseed = 101010
set.seed(iseed)

ate_before5_biased_prop=numeric()
ate_ipw5_biased_prop=numeric()
ate_aipw5_biased_prop=numeric()
aipw_se5_biased_prop=numeric()
aipw_ci_low5_biased_prop=numeric()
aipw_ci_high5_biased_prop=numeric()
sum_abs_error5_biased_prop = numeric()



for (i in seq_len(n_rep)) {
  
  data5 = list_data5[[i]]
  
  ate_before5_biased_prop[i] = mean(data5[which(data5$treat==0),]$dayslink) - mean(data5[which(data5$treat==1),]$dayslink)
  ate_ipw5_biased_prop[i] = ipw_biased(data5,formula)
  aipw_details = aipw_biased_prop(data5, outcome_fit = outcome_fit5[[i]], return_details = TRUE)
  ate_aipw5_biased_prop[i] = aipw_details$estimate
  aipw_se5_biased_prop[i] = aipw_details$se
  aipw_ci_low5_biased_prop[i] = aipw_details$ci_low
  aipw_ci_high5_biased_prop[i] = aipw_details$ci_high
  
  print(i)
}


ate5_biased_prop = c(mean(ate_before5_biased_prop),mean(ate_ipw5_biased_prop), mean(ate_aipw5_biased_prop))
var_ate5_biased_prop = c(var(ate_before5_biased_prop),var(ate_ipw5_biased_prop),var(ate_aipw5_biased_prop))


sum_abs_error5_biased_prop[1] = mean(abs(ate_before5_biased_prop-target5))
sum_abs_error5_biased_prop[2] = mean(abs(ate_ipw5_biased_prop-target5))
sum_abs_error5_biased_prop[3] = mean(abs(ate_aipw5_biased_prop-target5))



# Monte Carlo simulation intervals (diagnostic only)

as.numeric(quantile(ate_ipw2_biased_prop,c(0.025,0.975)))
as.numeric(quantile(ate_aipw2_biased_prop,c(0.025,0.975)))

as.numeric(quantile(ate_ipw3_biased_prop,c(0.025,0.975)))
as.numeric(quantile(ate_aipw3_biased_prop,c(0.025,0.975)))

as.numeric(quantile(ate_ipw4_biased_prop,c(0.025,0.975)))
as.numeric(quantile(ate_aipw4_biased_prop,c(0.025,0.975)))

as.numeric(quantile(ate_ipw5_biased_prop,c(0.025,0.975)))
as.numeric(quantile(ate_aipw5_biased_prop,c(0.025,0.975)))




####################################################################################



####################################################################################





ate_diff2 = abs(ate_before2-target2)
ate_diff3 = abs(ate_before3-target3)
ate_diff5 = abs(ate_before5-target5)


ate_plot2 = data.frame("x" = overlapping2, "y" = ate_diff2)
ate_plot3 = data.frame("x" = overlapping3, "y" = ate_diff3)
ate_plot5 = data.frame("x" = overlapping5, "y" = ate_diff5)

ggplot(data=ate_plot2, aes(x = x, y = y)) + geom_line(color ="blue") + geom_point(pch=16) + 
  
  xlab("non overlapping") + ylab("ATE difference") + ggtitle("Randomised new dataset")


ggplot(data=ate_plot3, aes(x = x, y = y)) + geom_line(color ="green") + geom_point(pch=16) + 
  
  xlab("non overlapping") + ylab("ATE difference") + ggtitle("Decoupling Mahalanobis")


ggplot(data=ate_plot5, aes(x = x, y = y)) + geom_line(color ="red") + geom_point(pch=16) + 
  
  xlab("non overlapping") + ylab("ATE difference") + ggtitle("Decoupling covariate")


xtable(ate)

xtable(sum_abs_error)


####################################################################################
# SAVE PAPER ANALYSIS OUTPUTS
####################################################################################

dir.create("results", showWarnings = FALSE, recursive = TRUE)

write.csv(ate, "results/HELP_table2_mean_estimates.csv", row.names = TRUE)
write.csv(sum_abs_error, "results/HELP_table2_mean_absolute_errors.csv", row.names = TRUE)

# The quantiles across the 1000 Monte Carlo estimates are simulation intervals,
# not confidence intervals for individual replicates.  Individual AIPW 95% CIs
# and their empirical coverage are calculated below from the AIPW influence function.

outcome_misspec_summary = data.frame(
  ITTE_method = c("Additional Sampling","Decoupling Mahalanobis","Perfect Knowledge","Decoupling Covariate"),
  Mean_target = c(mean(target2),mean(target3),mean(target4),mean(target5)),
  Unadjusted = c(mean(ate_before2_biased_outcome),mean(ate_before3_biased_outcome),mean(ate_before4_biased_outcome),mean(ate_before5_biased_outcome)),
  RB = c(mean(ate_reg2_biased_outcome),mean(ate_reg3_biased_outcome),mean(ate_reg4_biased_outcome),mean(ate_reg5_biased_outcome)),
  RB_MAE = c(sum_abs_error2_biased_outcome[2],sum_abs_error3_biased_outcome[2],sum_abs_error4_biased_outcome[2],sum_abs_error5_biased_outcome[2]),
  RB_simulation_interval_low = c(quantile(ate_reg2_biased_outcome,.025),quantile(ate_reg3_biased_outcome,.025),quantile(ate_reg4_biased_outcome,.025),quantile(ate_reg5_biased_outcome,.025)),
  RB_simulation_interval_high = c(quantile(ate_reg2_biased_outcome,.975),quantile(ate_reg3_biased_outcome,.975),quantile(ate_reg4_biased_outcome,.975),quantile(ate_reg5_biased_outcome,.975)),
  AIPW = c(mean(ate_aipw2_biased_outcome),mean(ate_aipw3_biased_outcome),mean(ate_aipw4_biased_outcome),mean(ate_aipw5_biased_outcome)),
  AIPW_MAE = c(sum_abs_error2_biased_outcome[3],sum_abs_error3_biased_outcome[3],sum_abs_error4_biased_outcome[3],sum_abs_error5_biased_outcome[3]),
  AIPW_simulation_interval_low = c(quantile(ate_aipw2_biased_outcome,.025),quantile(ate_aipw3_biased_outcome,.025),quantile(ate_aipw4_biased_outcome,.025),quantile(ate_aipw5_biased_outcome,.025)),
  AIPW_simulation_interval_high = c(quantile(ate_aipw2_biased_outcome,.975),quantile(ate_aipw3_biased_outcome,.975),quantile(ate_aipw4_biased_outcome,.975),quantile(ate_aipw5_biased_outcome,.975)),
  AIPW_mean_SE = c(mean(aipw_se2_biased_outcome),mean(aipw_se3_biased_outcome),mean(aipw_se4_biased_outcome),mean(aipw_se5_biased_outcome)),
  AIPW_mean_CI_width = c(mean(aipw_ci_high2_biased_outcome-aipw_ci_low2_biased_outcome),mean(aipw_ci_high3_biased_outcome-aipw_ci_low3_biased_outcome),mean(aipw_ci_high4_biased_outcome-aipw_ci_low4_biased_outcome),mean(aipw_ci_high5_biased_outcome-aipw_ci_low5_biased_outcome)),
  AIPW_coverage = c(mean(aipw_ci_low2_biased_outcome <= target2 & target2 <= aipw_ci_high2_biased_outcome),mean(aipw_ci_low3_biased_outcome <= target3 & target3 <= aipw_ci_high3_biased_outcome),mean(aipw_ci_low4_biased_outcome <= target4 & target4 <= aipw_ci_high4_biased_outcome),mean(aipw_ci_low5_biased_outcome <= target5 & target5 <= aipw_ci_high5_biased_outcome))
)

ps_misspec_summary = data.frame(
  ITTE_method = c("Additional Sampling","Decoupling Mahalanobis","Perfect Knowledge","Decoupling Covariate"),
  Mean_target = c(mean(target2),mean(target3),mean(target4),mean(target5)),
  Unadjusted = c(mean(ate_before2_biased_prop),mean(ate_before3_biased_prop),mean(ate_before4_biased_prop),mean(ate_before5_biased_prop)),
  IPW = c(mean(ate_ipw2_biased_prop),mean(ate_ipw3_biased_prop),mean(ate_ipw4_biased_prop),mean(ate_ipw5_biased_prop)),
  IPW_MAE = c(sum_abs_error2_biased_prop[2],sum_abs_error3_biased_prop[2],sum_abs_error4_biased_prop[2],sum_abs_error5_biased_prop[2]),
  IPW_simulation_interval_low = c(quantile(ate_ipw2_biased_prop,.025),quantile(ate_ipw3_biased_prop,.025),quantile(ate_ipw4_biased_prop,.025),quantile(ate_ipw5_biased_prop,.025)),
  IPW_simulation_interval_high = c(quantile(ate_ipw2_biased_prop,.975),quantile(ate_ipw3_biased_prop,.975),quantile(ate_ipw4_biased_prop,.975),quantile(ate_ipw5_biased_prop,.975)),
  AIPW = c(mean(ate_aipw2_biased_prop),mean(ate_aipw3_biased_prop),mean(ate_aipw4_biased_prop),mean(ate_aipw5_biased_prop)),
  AIPW_MAE = c(sum_abs_error2_biased_prop[3],sum_abs_error3_biased_prop[3],sum_abs_error4_biased_prop[3],sum_abs_error5_biased_prop[3]),
  AIPW_simulation_interval_low = c(quantile(ate_aipw2_biased_prop,.025),quantile(ate_aipw3_biased_prop,.025),quantile(ate_aipw4_biased_prop,.025),quantile(ate_aipw5_biased_prop,.025)),
  AIPW_simulation_interval_high = c(quantile(ate_aipw2_biased_prop,.975),quantile(ate_aipw3_biased_prop,.975),quantile(ate_aipw4_biased_prop,.975),quantile(ate_aipw5_biased_prop,.975)),
  AIPW_mean_SE = c(mean(aipw_se2_biased_prop),mean(aipw_se3_biased_prop),mean(aipw_se4_biased_prop),mean(aipw_se5_biased_prop)),
  AIPW_mean_CI_width = c(mean(aipw_ci_high2_biased_prop-aipw_ci_low2_biased_prop),mean(aipw_ci_high3_biased_prop-aipw_ci_low3_biased_prop),mean(aipw_ci_high4_biased_prop-aipw_ci_low4_biased_prop),mean(aipw_ci_high5_biased_prop-aipw_ci_low5_biased_prop)),
  AIPW_coverage = c(mean(aipw_ci_low2_biased_prop <= target2 & target2 <= aipw_ci_high2_biased_prop),mean(aipw_ci_low3_biased_prop <= target3 & target3 <= aipw_ci_high3_biased_prop),mean(aipw_ci_low4_biased_prop <= target4 & target4 <= aipw_ci_high4_biased_prop),mean(aipw_ci_low5_biased_prop <= target5 & target5 <= aipw_ci_high5_biased_prop))
)

standard_aipw_summary = data.frame(
  ITTE_method = c("Additional Sampling","Decoupling Mahalanobis","Perfect Knowledge","Decoupling Covariate"),
  Mean_target = c(mean(target2),mean(target3),mean(target4),mean(target5)),
  AIPW = c(mean(ate_aipw2),mean(ate_aipw3),mean(ate_aipw4),mean(ate_aipw5)),
  AIPW_MAE = c(sum_abs_error2[5],sum_abs_error3[5],sum_abs_error4[5],sum_abs_error5[5]),
  AIPW_mean_SE = c(mean(aipw_se2),mean(aipw_se3),mean(aipw_se4),mean(aipw_se5)),
  AIPW_mean_CI_width = c(mean(aipw_ci_high2-aipw_ci_low2),mean(aipw_ci_high3-aipw_ci_low3),mean(aipw_ci_high4-aipw_ci_low4),mean(aipw_ci_high5-aipw_ci_low5)),
  AIPW_coverage = c(mean(aipw_ci_low2 <= target2 & target2 <= aipw_ci_high2),mean(aipw_ci_low3 <= target3 & target3 <= aipw_ci_high3),mean(aipw_ci_low4 <= target4 & target4 <= aipw_ci_high4),mean(aipw_ci_low5 <= target5 & target5 <= aipw_ci_high5))
)

make_aipw_ci_rows <- function(scenario, specification, target, estimate, se, low, high) {
  data.frame(
    scenario = scenario,
    specification = specification,
    replicate = seq_along(target),
    target = target,
    estimate = estimate,
    se = se,
    ci_low = low,
    ci_high = high,
    covered = low <= target & target <= high
  )
}

aipw_ci_replicates = rbind(
  make_aipw_ci_rows("Additional Sampling","standard",target2,ate_aipw2,aipw_se2,aipw_ci_low2,aipw_ci_high2),
  make_aipw_ci_rows("Decoupling Mahalanobis","standard",target3,ate_aipw3,aipw_se3,aipw_ci_low3,aipw_ci_high3),
  make_aipw_ci_rows("Perfect Knowledge","standard",target4,ate_aipw4,aipw_se4,aipw_ci_low4,aipw_ci_high4),
  make_aipw_ci_rows("Decoupling Covariate","standard",target5,ate_aipw5,aipw_se5,aipw_ci_low5,aipw_ci_high5),
  make_aipw_ci_rows("Additional Sampling","outcome_misspecified",target2,ate_aipw2_biased_outcome,aipw_se2_biased_outcome,aipw_ci_low2_biased_outcome,aipw_ci_high2_biased_outcome),
  make_aipw_ci_rows("Decoupling Mahalanobis","outcome_misspecified",target3,ate_aipw3_biased_outcome,aipw_se3_biased_outcome,aipw_ci_low3_biased_outcome,aipw_ci_high3_biased_outcome),
  make_aipw_ci_rows("Perfect Knowledge","outcome_misspecified",target4,ate_aipw4_biased_outcome,aipw_se4_biased_outcome,aipw_ci_low4_biased_outcome,aipw_ci_high4_biased_outcome),
  make_aipw_ci_rows("Decoupling Covariate","outcome_misspecified",target5,ate_aipw5_biased_outcome,aipw_se5_biased_outcome,aipw_ci_low5_biased_outcome,aipw_ci_high5_biased_outcome),
  make_aipw_ci_rows("Additional Sampling","PS_misspecified",target2,ate_aipw2_biased_prop,aipw_se2_biased_prop,aipw_ci_low2_biased_prop,aipw_ci_high2_biased_prop),
  make_aipw_ci_rows("Decoupling Mahalanobis","PS_misspecified",target3,ate_aipw3_biased_prop,aipw_se3_biased_prop,aipw_ci_low3_biased_prop,aipw_ci_high3_biased_prop),
  make_aipw_ci_rows("Perfect Knowledge","PS_misspecified",target4,ate_aipw4_biased_prop,aipw_se4_biased_prop,aipw_ci_low4_biased_prop,aipw_ci_high4_biased_prop),
  make_aipw_ci_rows("Decoupling Covariate","PS_misspecified",target5,ate_aipw5_biased_prop,aipw_se5_biased_prop,aipw_ci_low5_biased_prop,aipw_ci_high5_biased_prop)
)

aipw_coverage_summary = aggregate(
  cbind(covered = as.numeric(aipw_ci_replicates$covered), ci_width = aipw_ci_replicates$ci_high-aipw_ci_replicates$ci_low, se = aipw_ci_replicates$se),
  by = list(scenario = aipw_ci_replicates$scenario, specification = aipw_ci_replicates$specification),
  FUN = mean
)
names(aipw_coverage_summary)[names(aipw_coverage_summary)=="covered"] = "coverage"
names(aipw_coverage_summary)[names(aipw_coverage_summary)=="ci_width"] = "mean_CI_width"
names(aipw_coverage_summary)[names(aipw_coverage_summary)=="se"] = "mean_SE"

write.csv(outcome_misspec_summary, "results/HELP_table3_misspecified_outcome.csv", row.names = FALSE)
write.csv(ps_misspec_summary, "results/HELP_table3_misspecified_PS.csv", row.names = FALSE)
write.csv(standard_aipw_summary, "results/HELP_table2_AIPW_CI_coverage.csv", row.names = FALSE)
write.csv(aipw_ci_replicates, "results/HELP_AIPW_CI_replicates.csv", row.names = FALSE)
write.csv(aipw_coverage_summary, "results/HELP_AIPW_CI_coverage_summary.csv", row.names = FALSE)

cat("\nAIPW 95% CI coverage (replicate-level influence-function intervals):\n")
print(aipw_coverage_summary, row.names = FALSE)

run_settings = data.frame(
  setting = c("n_rep","additional_q","additional_n_new","decoupling_mahalanobis_N","perfect_knowledge_gamma","perfect_knowledge_n","decoupling_covariate_N","correct_outcome_nuisance","PS_misspecified_omits"),
  value = c(as.character(n_rep),"0.90","60","150",as.character(perfect_gamma),"1000","90","lower-truncated Normal","female + daysdrink + daysanysub")
)
write.csv(run_settings, "results/HELP_run_settings.csv", row.names = FALSE)

replicate_targets = data.frame(
  replicate = seq_len(n_rep),
  Additional_Sampling = target2,
  Decoupling_Mahalanobis = target3,
  Perfect_Knowledge = target4,
  Decoupling_Covariate = target5
)
write.csv(replicate_targets, "results/HELP_replicate_targets.csv", row.names = FALSE)

target_summary = data.frame(
  scenario = c("Additional Sampling","Decoupling Mahalanobis","Perfect Knowledge","Decoupling Covariate"),
  mean_target = c(mean(target2),mean(target3),mean(target4),mean(target5)),
  sd_target = c(sd(target2),sd(target3),sd(target4),sd(target5))
)
write.csv(target_summary, "results/HELP_target_summary.csv", row.names = FALSE)
cat("\nReplicate-specific target summary:\n")
print(target_summary, row.names = FALSE)

save(
  ate, sum_abs_error,
  target2, target3, target4, target5, replicate_targets, target_summary,
  outcome_misspec_summary, ps_misspec_summary, standard_aipw_summary, aipw_ci_replicates, aipw_coverage_summary,
  overlapping2, overlapping3, overlapping4, overlapping5,
  ate_before2, ate_before3, ate_before4, ate_before5,
  att_tte2, att_tte3, att_tte4, att_tte5,
  ate_ipw2, ate_ipw3, ate_ipw4, ate_ipw5,
  ate_reg2, ate_reg3, ate_reg4, ate_reg5,
  ate_aipw2, ate_aipw3, ate_aipw4, ate_aipw5,
  aipw_se2, aipw_se3, aipw_se4, aipw_se5,
  aipw_ci_low2, aipw_ci_low3, aipw_ci_low4, aipw_ci_low5,
  aipw_ci_high2, aipw_ci_high3, aipw_ci_high4, aipw_ci_high5,
  ate_reg2_biased_outcome, ate_reg3_biased_outcome, ate_reg4_biased_outcome, ate_reg5_biased_outcome,
  ate_aipw2_biased_outcome, ate_aipw3_biased_outcome, ate_aipw4_biased_outcome, ate_aipw5_biased_outcome,
  aipw_se2_biased_outcome, aipw_se3_biased_outcome, aipw_se4_biased_outcome, aipw_se5_biased_outcome,
  aipw_ci_low2_biased_outcome, aipw_ci_low3_biased_outcome, aipw_ci_low4_biased_outcome, aipw_ci_low5_biased_outcome,
  aipw_ci_high2_biased_outcome, aipw_ci_high3_biased_outcome, aipw_ci_high4_biased_outcome, aipw_ci_high5_biased_outcome,
  ate_ipw2_biased_prop, ate_ipw3_biased_prop, ate_ipw4_biased_prop, ate_ipw5_biased_prop,
  ate_aipw2_biased_prop, ate_aipw3_biased_prop, ate_aipw4_biased_prop, ate_aipw5_biased_prop,
  aipw_se2_biased_prop, aipw_se3_biased_prop, aipw_se4_biased_prop, aipw_se5_biased_prop,
  aipw_ci_low2_biased_prop, aipw_ci_low3_biased_prop, aipw_ci_low4_biased_prop, aipw_ci_low5_biased_prop,
  aipw_ci_high2_biased_prop, aipw_ci_high3_biased_prop, aipw_ci_high4_biased_prop, aipw_ci_high5_biased_prop,
  file = "results/HELP_paper_analysis_results.RData"
)

cat("\nPaper-analysis outputs saved in ./results/\n")


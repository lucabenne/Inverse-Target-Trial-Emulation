dir.create("results", showWarnings=FALSE, recursive=TRUE)

library(TH.data)
library(survival)
library(survHE)
library(R2jags)
library(truncnorm)

cat("\n================ PACKAGE VERSIONS ================\n")
cat("R:", R.version.string, "\n")
cat("survHE:", as.character(packageVersion("survHE")), "\n")
cat("survHEhmc:", as.character(packageVersion("survHEhmc")), "\n")
cat("rstan:", as.character(packageVersion("rstan")), "\n")

source("../src/00_GBSG2_helpers.R")
source("../src/03_fit_GBSG2_covariate_DGP.R")
source("../src/05_generate_GBSG2_synthetic_trial.R")
source("../src/10_estimate_weibull_survival_model.R")

# -----------------------------------------------------------------------------
# Source RCT: same recoding as the original GBSG2 scripts
# -----------------------------------------------------------------------------
data("GBSG2")
datt = GBSG2

datt$tgrade=as.numeric(datt$tgrade)
datt$horTh=as.integer(as.integer(datt$horTh)-1)
datt$menostat=as.integer(as.integer(datt$menostat)-1)
colnames(datt)[10]="event"
colnames(datt)[1]="treat"

# -----------------------------------------------------------------------------
# Source-based standardisation for numerical stability in survHEhmc only.
# -----------------------------------------------------------------------------
scaling = get_GBSG2_survival_scaling(datt)
datt_model = apply_GBSG2_survival_scaling(datt,scaling)

formula_full = Surv(time,event) ~ as.factor(treat) +
  age_s + menostat_s + tsize_s + tgrade_s + pnodes_s + progrec_s + estrec_s

distrib = "weibull"
iseed = 100001
set.seed(iseed)

cat("\n>>> SOURCE WEIBULL EVENT + CENSORING FIT START\n")
weib = estimate_params(datt_model,formula_full,distrib,msd=0,seed=iseed)
cat("\n>>> SOURCE WEIBULL EVENT + CENSORING FIT DONE\n")

cat("\n================ SOURCE WEIBULL PARAMETERS ================\n")
cat("Event shape alpha     =",weib$alpha_event,"\n")
cat("Censoring shape alpha =",weib$alpha_cens,"\n")
cat("Source adjusted MSD   =",weib$mean_surv_diff,"\n")

# Reset seed before posterior-predictive synthetic-RCT generation.
set.seed(iseed)
cat("\n>>> SYNTHETIC RCT GENERATION START\n")
synthetic_rct = generate_trial(
  datt,
  weib$beta_event,
  weib$beta_cens,
  weib$alpha_event,
  weib$alpha_cens,
  scaling
)
cat("\n>>> SYNTHETIC RCT GENERATION DONE\n")

# -----------------------------------------------------------------------------
# 1. Marginal covariate distributions
# -----------------------------------------------------------------------------
covs = c("age","menostat","tsize","tgrade","pnodes","progrec","estrec")

marginal_table = do.call(rbind,lapply(covs,function(v){
  x = datt[[v]]
  y = synthetic_rct[[v]]
  pooled_sd = sqrt((var(x)+var(y))/2)
  std_diff = ifelse(pooled_sd>0,(mean(y)-mean(x))/pooled_sd,NA)
  data.frame(
    variable=v,
    source_mean=mean(x),
    synthetic_mean=mean(y),
    source_sd=sd(x),
    synthetic_sd=sd(y),
    standardised_difference=std_diff
  )
}))

cat("\n================ MARGINAL COVARIATES ================\n")
print(marginal_table,row.names=FALSE,digits=5)
cat("\nMean |standardised difference| =",mean(abs(marginal_table$standardised_difference),na.rm=TRUE),"\n")
cat("Max  |standardised difference| =",max(abs(marginal_table$standardised_difference),na.rm=TRUE),"\n")

marginal_table$sd_ratio_synthetic_source = marginal_table$synthetic_sd / marginal_table$source_sd
cat("\n================ SD RATIOS ================\n")
print(marginal_table[,c("variable","source_sd","synthetic_sd","sd_ratio_synthetic_source")],
      row.names=FALSE,digits=5)

cat("\n================ TGRADE SUPPORT ================\n")
cat("Source counts:\n")
print(table(datt$tgrade))
cat("Synthetic counts:\n")
print(table(synthetic_rct$tgrade))
if(any(!synthetic_rct$tgrade %in% c(1,2,3))){
  stop("Synthetic tgrade contains values outside {1,2,3}.")
}

shape_table = do.call(rbind,lapply(covs,function(v){
  x=datt[[v]]; y=synthetic_rct[[v]]
  data.frame(
    variable=v,
    source_median=median(x), synthetic_median=median(y),
    source_q25=as.numeric(quantile(x,.25)), synthetic_q25=as.numeric(quantile(y,.25)),
    source_q75=as.numeric(quantile(x,.75)), synthetic_q75=as.numeric(quantile(y,.75)),
    source_q975=as.numeric(quantile(x,.975)), synthetic_q975=as.numeric(quantile(y,.975))
  )
}))
cat("\n================ MARGINAL SHAPE CHECK ================\n")
print(shape_table,row.names=FALSE,digits=5)

# -----------------------------------------------------------------------------
# 2. Dependence structure
# -----------------------------------------------------------------------------
cor_source = cor(datt[,covs],use="pairwise.complete.obs")
cor_synthetic = cor(synthetic_rct[,covs],use="pairwise.complete.obs")
cor_diff = cor_synthetic-cor_source

upper = upper.tri(cor_diff)
cor_long = data.frame(
  var1=row.names(cor_diff)[row(cor_diff)[upper]],
  var2=colnames(cor_diff)[col(cor_diff)[upper]],
  source=cor_source[upper],
  synthetic=cor_synthetic[upper],
  difference=cor_diff[upper]
)
cor_long$abs_difference=abs(cor_long$difference)
cor_long=cor_long[order(-cor_long$abs_difference),]

cat("\n================ CORRELATIONS ================\n")
cat("Mean |correlation difference| =",mean(cor_long$abs_difference),"\n")
cat("Max  |correlation difference| =",max(cor_long$abs_difference),"\n")
cat("\nTop 10 absolute correlation discrepancies:\n")
print(head(cor_long,10),row.names=FALSE,digits=5)

# -----------------------------------------------------------------------------
# 3. Outcome distribution
# -----------------------------------------------------------------------------
outcome_table = data.frame(
  Dataset=c("Source RCT","Synthetic RCT"),
  Mean_followup=c(mean(datt$time),mean(synthetic_rct$time)),
  SD_followup=c(sd(datt$time),sd(synthetic_rct$time)),
  Median_followup=c(median(datt$time),median(synthetic_rct$time)),
  Q25=c(as.numeric(quantile(datt$time,.25)),as.numeric(quantile(synthetic_rct$time,.25))),
  Q75=c(as.numeric(quantile(datt$time,.75)),as.numeric(quantile(synthetic_rct$time,.75))),
  Q975=c(as.numeric(quantile(datt$time,.975)),as.numeric(quantile(synthetic_rct$time,.975))),
  Event_proportion=c(mean(datt$event),mean(synthetic_rct$event)),
  Censoring_proportion=c(mean(1-datt$event),mean(1-synthetic_rct$event))
)

cat("\n================ OUTCOME DISTRIBUTION ================\n")
print(outcome_table,row.names=FALSE,digits=5)

# -----------------------------------------------------------------------------
# 4. Treatment allocation -- design check only
# -----------------------------------------------------------------------------
treatment_table = data.frame(
  Dataset=c("Source RCT","Synthetic RCT"),
  Treated_proportion=c(mean(datt$treat),mean(synthetic_rct$treat))
)
cat("\n================ TREATMENT ALLOCATION ================\n")
print(treatment_table,row.names=FALSE,digits=5)

# -----------------------------------------------------------------------------
# 5. Adjusted treatment effect under the SAME Weibull AFT model
# -----------------------------------------------------------------------------
source_effect = weib$mean_surv_diff
synthetic_model = apply_GBSG2_survival_scaling(synthetic_rct,scaling)

cat("\n>>> SYNTHETIC WEIBULL TREATMENT-EFFECT FIT START\n")
synthetic_fit = estimate_params(
  synthetic_model,formula_full,distrib,msd=1,seed=iseed+100L
)
synthetic_effect = synthetic_fit$mean_surv_diff
cat("\n>>> SYNTHETIC WEIBULL TREATMENT-EFFECT FIT DONE\n")

treatment_effect_validation = data.frame(
  Dataset=c("Source RCT","Synthetic RCT"),
  Adjusted_mean_survival_difference=c(source_effect,synthetic_effect),
  Event_shape_alpha=c(weib$alpha_event,synthetic_fit$alpha_event)
)

cat("\n================ ADJUSTED TREATMENT EFFECT ================\n")
print(treatment_effect_validation,row.names=FALSE,digits=5)

# -----------------------------------------------------------------------------
# 6. SOURCE-ONLY Weibull posterior-predictive replication check.
# Same source patients/covariates/treatment; only event/censor times regenerated.
# This isolates adequacy of the chosen survival/censoring family.
# -----------------------------------------------------------------------------
set.seed(20261004)
B_ppc = 1000
X_source = model.matrix(
  ~ as.factor(treat) + age_s + menostat_s + tsize_s + tgrade_s +
    pnodes_s + progrec_s + estrec_s,
  data=datt_model
)

scale_event_source = exp(as.numeric(X_source %*% weib$beta_event))
scale_cens_source = exp(as.numeric(X_source %*% weib$beta_cens))

ppc = matrix(NA_real_,B_ppc,7)
colnames(ppc)=c("mean_followup","sd_followup","median_followup","q25","q75","q975","event_proportion")

for(b in seq_len(B_ppc)){
  te = rweibull(nrow(datt),shape=weib$alpha_event,scale=scale_event_source)
  tc = rweibull(nrow(datt),shape=weib$alpha_cens,scale=scale_cens_source)
  tobs = pmin(te,tc)
  eobs = as.numeric(te<tc)
  ppc[b,] = c(
    mean(tobs),sd(tobs),median(tobs),
    as.numeric(quantile(tobs,.25)),
    as.numeric(quantile(tobs,.75)),
    as.numeric(quantile(tobs,.975)),
    mean(eobs)
  )
}
ppc = as.data.frame(ppc)

observed_ppc = c(
  mean_followup=mean(datt$time),
  sd_followup=sd(datt$time),
  median_followup=median(datt$time),
  q25=as.numeric(quantile(datt$time,.25)),
  q75=as.numeric(quantile(datt$time,.75)),
  q975=as.numeric(quantile(datt$time,.975)),
  event_proportion=mean(datt$event)
)

ppc_summary = do.call(rbind,lapply(names(observed_ppc),function(v){
  sims=ppc[[v]]; obs=observed_ppc[v]
  data.frame(
    Quantity=v,
    Source_observed=as.numeric(obs),
    Replicated_mean=mean(sims),
    Replicated_q025=as.numeric(quantile(sims,.025)),
    Replicated_q975=as.numeric(quantile(sims,.975)),
    Inside_95pct_PPC=(obs>=quantile(sims,.025) & obs<=quantile(sims,.975))
  )
}))
rownames(ppc_summary)=NULL

cat("\n================ SOURCE-ONLY WEIBULL PPC ================\n")
print(ppc_summary,row.names=FALSE,digits=5)

cat("\n================ COMPACT CHECK ================\n")
cat("Mean |standardised difference|:",mean(abs(marginal_table$standardised_difference),na.rm=TRUE),"\n")
cat("Max  |standardised difference|:",max(abs(marginal_table$standardised_difference),na.rm=TRUE),"\n")
cat("Mean |correlation difference|:",mean(cor_long$abs_difference),"\n")
cat("Max  |correlation difference|:",max(cor_long$abs_difference),"\n")
cat("Source event proportion:",mean(datt$event),"\n")
cat("Synthetic event proportion:",mean(synthetic_rct$event),"\n")
cat("Source adjusted Weibull MSD:",source_effect,"\n")
cat("Synthetic adjusted Weibull MSD:",synthetic_effect,"\n")
cat("Source event alpha:",weib$alpha_event,"\n")
cat("Synthetic event alpha:",synthetic_fit$alpha_event,"\n")

# Save validation outputs for the revision record.
write.csv(marginal_table,"results/GBSG2_validation_covariate_marginals_weibull.csv",row.names=FALSE)
write.csv(shape_table,"results/GBSG2_validation_covariate_shapes_weibull.csv",row.names=FALSE)
write.csv(cor_long,"results/GBSG2_validation_correlations_weibull.csv",row.names=FALSE)
write.csv(outcome_table,"results/GBSG2_validation_outcome_weibull.csv",row.names=FALSE)
write.csv(treatment_effect_validation,"results/GBSG2_validation_treatment_effect_weibull.csv",row.names=FALSE)
write.csv(ppc_summary,"results/GBSG2_validation_source_PPC_weibull.csv",row.names=FALSE)
save(
  synthetic_rct,marginal_table,shape_table,cor_source,cor_synthetic,cor_long,
  outcome_table,treatment_effect_validation,ppc_summary,weib,scaling,
  file="results/GBSG2_SYNTHETIC_RCT_VALIDATION_WEIBULL.RData"
)

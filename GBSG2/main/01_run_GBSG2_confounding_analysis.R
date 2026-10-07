############################################################
# 01_run_GBSG2_confounding_analysis.R
#
# Revised GBSG2 confounding experiment.
#
# Key Monte Carlo design:
#   1) Fit the source Weibull event/censoring DGP ONCE.
#   2) Fit the covariate JAGS DGP (BG) ONCE.
#   3) Use fixed severity levels N = 20,40,60,80,100.
#   4) For every Monte Carlo replicate, draw ONE joint posterior
#      parameter vector from the already-fitted BG.
#   5) Decoupling: generate a fresh randomized synthetic RCT from that
#      posterior draw, then deterministically apply each N level.
#   6) Additional Sampling: keep the original source-RCT-based mechanism,
#      but use one BG posterior draw per replicate for the local conditional
#      generation and new stochastic additions at each replicate.
#
# NO JAGS or source-Weibull refit occurs inside the replication loops.
############################################################

library(TH.data)
library(survival)
library(R2jags)
library(truncnorm)
library(ggplot2)

source("../src/00_GBSG2_helpers.R")
source("../src/03_fit_GBSG2_covariate_DGP.R")
source("../src/05_generate_GBSG2_synthetic_trial.R")
source("../src/06_generate_additional_sampling_survival.R")
source("../src/07_generate_decoupling_covariate_survival.R")
source("../src/08_measure_mahalanobis_imbalance_survival.R")
source("../src/09_pairwise_mahalanobis_survival.R")
source("../src/10_estimate_weibull_survival_model.R")

############################################################
# USER-CONTROLLED SETTINGS
############################################################

iseed <- 100001L
N_levels <- c(20L,40L,60L,80L,100L)
n_rep <- 30L
q_additional <- 0.80
x_decoupling <- "pnodes"
nsim <- 1000L

# FALSE avoids saving 300 large data frames.  Set TRUE only for debugging.
store_datasets <- FALSE

if(any(diff(N_levels)<=0L)) stop("N_levels must be strictly increasing.")
if(max(N_levels)>=nsim/2) stop("N_levels are too large relative to nsim.")

############################################################
# SOURCE DATA + FIT DGP ONCE
############################################################

set.seed(iseed)
datt <- prepare_GBSG2_data()
scaling <- get_GBSG2_survival_scaling(datt)
datt_model <- apply_GBSG2_survival_scaling(datt,scaling)

formula_full <- GBSG2_full_survival_formula()
formula_naive <- Surv(time,event) ~ as.factor(treat)
distrib <- "weibull"

cat("\n============================================================\n")
cat("GBSG2 CONFOUNDING ANALYSIS -- POSTERIOR REPLICATES\n")
cat("============================================================\n")
cat("N levels:",paste(N_levels,collapse=", "),"\n")
cat("Monte Carlo replicates:",n_rep,"\n")
cat("Synthetic RCT size:",nsim,"\n\n")

cat(">>> FIT SOURCE WEIBULL EVENT/CENSORING DGP ONCE\n")
weib <- estimate_params(
  datt_model,formula_full,distrib,msd=0,seed=iseed
)

cat("\n>>> FIT COVARIATE JAGS DGP (BG) ONCE\n")
BG <- generate_BG(datt)
n_BG_draws <- GBSG2_BG_n_draws(BG)
cat("Available joint BG posterior draws:",n_BG_draws,"\n")

# Source naive treatment-only reference for Additional Sampling.
cat("\n>>> FIT SOURCE NAIVE WEIBULL MSD ONCE\n")
source_naive_msd <- estimate_params(
  datt,formula_naive,distrib,msd=1,seed=iseed+10L
)$mean_surv_diff

cat("Source adjusted Weibull MSD =",round(weib$mean_surv_diff,3),"\n")
cat("Source naive Weibull MSD    =",round(source_naive_msd,3),"\n")

############################################################
# FIXED OBJECTS USED BY ADDITIONAL SAMPLING
############################################################

# Same source-covariate Mahalanobis definition as the original experiment.
index_additional <- c(1,2,3,9,10)
Ma <- pairwise_mahalanobis(datt,index_additional)
source_overlap <- overlapping_surv(datt,index_additional)

############################################################
# PRE-DRAW ONE JOINT BG INDEX PER REPLICATE / METHOD
#
# BG IS NOT REFITTED.  These are just indices into the already
# available posterior sample.
############################################################

set.seed(iseed+500L)
BG_draw_additional <- sample.int(n_BG_draws,n_rep,replace=TRUE)
BG_draw_decoupling <- sample.int(n_BG_draws,n_rep,replace=TRUE)

############################################################
# STORAGE
############################################################

K <- length(N_levels)

additional_overlap <- matrix(NA_real_,nrow=n_rep,ncol=K)
additional_msd <- matrix(NA_real_,nrow=n_rep,ncol=K)
additional_abs_change <- matrix(NA_real_,nrow=n_rep,ncol=K)
additional_n_final <- matrix(NA_integer_,nrow=n_rep,ncol=K)

decoupling_baseline_overlap <- rep(NA_real_,n_rep)
decoupling_baseline_msd <- rep(NA_real_,n_rep)
decoupling_overlap <- matrix(NA_real_,nrow=n_rep,ncol=K)
decoupling_msd <- matrix(NA_real_,nrow=n_rep,ncol=K)
decoupling_abs_change <- matrix(NA_real_,nrow=n_rep,ncol=K)
decoupling_n_final <- matrix(NA_integer_,nrow=n_rep,ncol=K)

colnames(additional_overlap) <- colnames(additional_msd) <-
  colnames(additional_abs_change) <- colnames(additional_n_final) <- paste0("N",N_levels)
colnames(decoupling_overlap) <- colnames(decoupling_msd) <-
  colnames(decoupling_abs_change) <- colnames(decoupling_n_final) <- paste0("N",N_levels)

ls_newmahal <- if(store_datasets) vector("list",n_rep) else NULL
ls_decoupling <- if(store_datasets) vector("list",n_rep) else NULL

############################################################
# 1. ADDITIONAL SAMPLING
#
# N levels are cumulative.  Example: for N=20,40,... we append
# 20 new observations at each step, but use a fresh Bernoulli
# sequence and a fresh BG posterior draw in every replicate.
############################################################

cat("\n============================================================\n")
cat("ADDITIONAL SAMPLING\n")
cat("============================================================\n")

N_increment <- diff(c(0L,N_levels))

for(j in seq_len(n_rep)){

  set.seed(iseed + 10000L + j)
  bern_j <- rbinom(max(N_levels),1,0.5)
  data_add <- datt
  previous_N <- 0L

  if(store_datasets) ls_newmahal[[j]] <- vector("list",K)

  cat("Additional Sampling replicate",j,"/",n_rep,
      " | BG draw",BG_draw_additional[j],"\n")

  for(k in seq_len(K)){

    n_add <- N_increment[k]
    idx_bern <- seq.int(previous_N+1L,N_levels[k])

    data_add <- newdatamahal_distr_surv2(
      data=datt,
      bern=bern_j[idx_bern],
      n=n_add,
      q=q_additional,
      beta_weib=weib$beta_event,
      beta_weib_cens=weib$beta_cens,
      alpha_event=weib$alpha_event,
      alpha_cens=weib$alpha_cens,
      scaling=scaling,
      Ma=Ma,
      BG=BG,
      data_prec=data_add,
      BG_draw_index=BG_draw_additional[j]
    )

    additional_overlap[j,k] <- overlapping_surv(data_add,index_additional)

    additional_msd[j,k] <- estimate_params(
      data_add,formula_naive,distrib,msd=1,
      seed=iseed + 100000L + j*100L + k
    )$mean_surv_diff

    additional_abs_change[j,k] <- abs(additional_msd[j,k]-source_naive_msd)
    additional_n_final[j,k] <- nrow(data_add)

    if(store_datasets) ls_newmahal[[j]][[k]] <- data_add
    previous_N <- N_levels[k]
  }
}

############################################################
# 2. DECOUPLING COVARIATE
#
# Every replicate gets a NEW synthetic RCT generated from one
# posterior BG draw.  N stays fixed at the pre-specified levels.
# The decoupling removal rule remains deterministic conditional
# on that synthetic RCT.
############################################################

cat("\n============================================================\n")
cat("DECOUPLING COVARIATE\n")
cat("============================================================\n")

index_decoupling <- c(1,9,10)

for(j in seq_len(n_rep)){

  cat("Decoupling replicate",j,"/",n_rep,
      " | BG draw",BG_draw_decoupling[j],"\n")

  set.seed(iseed + 20000L + j)

  sim_data <- generate_trial(
    data=datt,
    beta_weib=weib$beta_event,
    beta_weib_cens=weib$beta_cens,
    alpha_event=weib$alpha_event,
    alpha_cens=weib$alpha_cens,
    scaling=scaling,
    BG=BG,
    nsim=nsim,
    posterior_mode="dataset",
    posterior_draw_index=BG_draw_decoupling[j]
  )

  n0 <- sum(sim_data$treat==0)
  n1 <- sum(sim_data$treat==1)
  if(min(n0,n1)<=max(N_levels)){
    stop("Replicate ",j," has too few subjects in one treatment group for max N=",
         max(N_levels),". n0=",n0,", n1=",n1)
  }

  decoupling_baseline_overlap[j] <- overlapping_surv(sim_data,index_decoupling)
  decoupling_baseline_msd[j] <- estimate_params(
    sim_data,formula_naive,distrib,msd=1,
    seed=iseed + 200000L + j*100L
  )$mean_surv_diff

  if(store_datasets) ls_decoupling[[j]] <- vector("list",K)

  for(k in seq_len(K)){

    data_decoupling <- decoupling_covariate_surv1(
      N=N_levels[k],
      p=0.5,
      x=x_decoupling,
      beta_weib=weib$beta_event,
      beta_weib_cens=weib$beta_cens,
      sim_data=sim_data
    )

    decoupling_overlap[j,k] <- overlapping_surv(data_decoupling,index_decoupling)

    decoupling_msd[j,k] <- estimate_params(
      data_decoupling,formula_naive,distrib,msd=1,
      seed=iseed + 300000L + j*100L + k
    )$mean_surv_diff

    # Pair every biased dataset with ITS OWN pre-decoupling synthetic RCT.
    decoupling_abs_change[j,k] <-
      abs(decoupling_msd[j,k]-decoupling_baseline_msd[j])

    decoupling_n_final[j,k] <- nrow(data_decoupling)

    if(store_datasets) ls_decoupling[[j]][[k]] <- data_decoupling
  }
}

############################################################
# LONG-FORM REPLICATION OUTPUT
############################################################

make_long <- function(method,overlap,msd,abs_change,n_final,bg_draw,baseline_msd){
  out <- vector("list",n_rep*K)
  z <- 1L
  for(j in seq_len(n_rep)){
    for(k in seq_len(K)){
      out[[z]] <- data.frame(
        method=method,
        replicate=j,
        N=N_levels[k],
        BG_draw_index=bg_draw[j],
        n_final=n_final[j,k],
        baseline_MSD=if(length(baseline_msd)==1L) baseline_msd else baseline_msd[j],
        Mahalanobis=overlap[j,k],
        observational_MSD=msd[j,k],
        absolute_MSD_change=abs_change[j,k]
      )
      z <- z+1L
    }
  }
  do.call(rbind,out)
}

additional_long <- make_long(
  "Additional Sampling",additional_overlap,additional_msd,
  additional_abs_change,additional_n_final,BG_draw_additional,source_naive_msd
)

decoupling_long <- make_long(
  "Decoupling covariate",decoupling_overlap,decoupling_msd,
  decoupling_abs_change,decoupling_n_final,BG_draw_decoupling,decoupling_baseline_msd
)

confounding_replications <- rbind(additional_long,decoupling_long)

############################################################
# SUMMARY
############################################################

mean_col <- function(x) colMeans(x,na.rm=TRUE)
sd_col <- function(x) apply(x,2,sd,na.rm=TRUE)
q025_col <- function(x) apply(x,2,quantile,probs=0.025,na.rm=TRUE)
q975_col <- function(x) apply(x,2,quantile,probs=0.975,na.rm=TRUE)

confounding_summary <- data.frame(
  N=N_levels,
  additional_total_added=N_levels,
  decoupling_total_removed=2L*N_levels,

  additional_sampling_mahalanobis=mean_col(additional_overlap),
  additional_sampling_mahalanobis_sd=sd_col(additional_overlap),
  additional_sampling_mean_abs_MSD_change=mean_col(additional_abs_change),
  additional_sampling_abs_MSD_change_sd=sd_col(additional_abs_change),
  additional_sampling_abs_MSD_change_q025=q025_col(additional_abs_change),
  additional_sampling_abs_MSD_change_q975=q975_col(additional_abs_change),

  decoupling_mahalanobis=mean_col(decoupling_overlap),
  decoupling_mahalanobis_sd=sd_col(decoupling_overlap),
  decoupling_mean_abs_MSD_change=mean_col(decoupling_abs_change),
  decoupling_abs_MSD_change_sd=sd_col(decoupling_abs_change),
  decoupling_abs_MSD_change_q025=q025_col(decoupling_abs_change),
  decoupling_abs_MSD_change_q975=q975_col(decoupling_abs_change)
)

cat("\n============================================================\n")
cat("CONFOUNDING SUMMARY\n")
cat("============================================================\n")
print(confounding_summary,row.names=FALSE,digits=5)

cat("\nDecoupling baseline synthetic-RCT naive MSD:\n")
cat("  mean =",round(mean(decoupling_baseline_msd),3),"\n")
cat("  SD   =",round(sd(decoupling_baseline_msd),3),"\n")
cat("  range=",paste(round(range(decoupling_baseline_msd),3),collapse=" to "),"\n")

cat("\nDecoupling baseline Mahalanobis:\n")
cat("  mean =",round(mean(decoupling_baseline_overlap),4),"\n")
cat("  SD   =",round(sd(decoupling_baseline_overlap),4),"\n")

############################################################
# SAVE TABLES
############################################################

write.csv(
  confounding_summary,
  "GBSG2_confounding_summary_weibull_posterior_replicates.csv",
  row.names=FALSE
)

write.csv(
  confounding_replications,
  "GBSG2_confounding_replications_weibull.csv",
  row.names=FALSE
)

############################################################
# PLOTS: MEAN SEVERITY VS MEAN ABSOLUTE MSD CHANGE
############################################################

ate_plot_additional <- data.frame(
  x=mean_col(additional_overlap),
  y=mean_col(additional_abs_change)
)

ate_plot_decoupling <- data.frame(
  x=mean_col(decoupling_overlap),
  y=mean_col(decoupling_abs_change)
)

plot_additional <- ggplot(ate_plot_additional,aes(x=x,y=y))+
  geom_line()+geom_point()+
  xlab("Mean non-overlap (Mahalanobis distance)")+
  ylab("Mean absolute change in Weibull MSD")+
  ggtitle("Additional Sampling")

plot_decoupling <- ggplot(ate_plot_decoupling,aes(x=x,y=y))+
  geom_line()+geom_point()+
  xlab("Mean non-overlap (Mahalanobis distance)")+
  ylab("Mean absolute change in Weibull MSD")+
  ggtitle("Decoupling covariate")

print(plot_additional)
print(plot_decoupling)

############################################################
# SAVE R OBJECTS
############################################################

save(
  weib,BG,
  N_levels,n_rep,q_additional,x_decoupling,nsim,
  source_naive_msd,source_overlap,
  BG_draw_additional,BG_draw_decoupling,
  additional_overlap,additional_msd,additional_abs_change,additional_n_final,
  decoupling_baseline_overlap,decoupling_baseline_msd,
  decoupling_overlap,decoupling_msd,decoupling_abs_change,decoupling_n_final,
  confounding_summary,confounding_replications,
  ls_newmahal,ls_decoupling,
  file="GBSG2_confounding_results_weibull_posterior_replicates.RData"
)

cat("\n============================================================\n")
cat("DONE\n")
cat("BG was fitted once; source Weibull was fitted once.\n")
cat("Every decoupling replicate used a fresh synthetic RCT from one BG posterior draw.\n")
cat("N levels were fixed at:",paste(N_levels,collapse=", "),"\n")
cat("============================================================\n")

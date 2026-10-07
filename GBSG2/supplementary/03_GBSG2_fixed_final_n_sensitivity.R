###############################################################################
# 03_GBSG2_fixed_final_n_sensitivity.R
#
# Reviewer 3 sensitivity: vary the strength of the bias-generating mechanism
# while holding FINAL sample size fixed at n=1000.
#
# Additional Sampling:
#   create a randomized synthetic baseline of size 1000-N, then append N
#   source-anchored Additional-Sampling subjects -> final n=1000.
#
# Decoupling Covariate:
#   create a randomized synthetic baseline of size 1000+2N, then remove N
#   controls and N treated subjects -> final n=1000.
#
# Thus the final sample size is constant by construction; changes in SMD /
# PS overlap cannot be attributed to changing final n.
###############################################################################

if (!exists("weib", inherits = FALSE) || !exists("BG", inherits = FALSE)) {
  
}
source("99_GBSG2_diagnostic_helpers.R")

N_levels <- c(20L,40L,60L,80L,100L)
n_rep <- 30L
final_n <- 1000L
q_additional <- 0.80
x_decoupling <- "pnodes"

index_additional <- c(1,2,3,9,10)
Ma_source <- pairwise_mahalanobis(datt, index_additional)

set.seed(130501L)
BG_draws <- sample.int(n_BG_draws, n_rep, replace = TRUE)

rep_rows <- list()
zz <- 1L

for (j in seq_len(n_rep)) {
  draw_j <- BG_draws[j]

  for (k in seq_along(N_levels)) {
    Nk <- N_levels[k]

    # -------------------------
    # Additional Sampling
    # -------------------------
    set.seed(131000L + j*100L + k)
    baseline_add <- generate_trial(
      data = datt,
      beta_weib = weib$beta_event,
      beta_weib_cens = weib$beta_cens,
      alpha_event = weib$alpha_event,
      alpha_cens = weib$alpha_cens,
      scaling = scaling,
      BG = BG,
      nsim = final_n - Nk,
      posterior_mode = "dataset",
      posterior_draw_index = draw_j
    )

    bern <- rbinom(Nk, 1, 0.5)

    d_add <- newdatamahal_distr_surv2(
      data = datt,
      bern = bern,
      n = Nk,
      q = q_additional,
      beta_weib = weib$beta_event,
      beta_weib_cens = weib$beta_cens,
      alpha_event = weib$alpha_event,
      alpha_cens = weib$alpha_cens,
      scaling = scaling,
      Ma = Ma_source,
      BG = BG,
      data_prec = baseline_add,
      BG_draw_index = draw_j
    )
    if (nrow(d_add) != final_n) stop("Additional fixed-n construction failed.")

    one <- GBSG2_diag_one(d_add, source_data = datt, mahal_index = index_additional)
    rr <- as.data.frame(as.list(one))
    rr$method <- "Additional Sampling"
    rr$N <- Nk
    rr$replicate <- j
    rr$final_n_target <- final_n
    rep_rows[[zz]] <- rr
    zz <- zz + 1L

    # -------------------------
    # Decoupling Covariate
    # -------------------------
    set.seed(132000L + j*100L + k)
    baseline_dec <- generate_trial(
      data = datt,
      beta_weib = weib$beta_event,
      beta_weib_cens = weib$beta_cens,
      alpha_event = weib$alpha_event,
      alpha_cens = weib$alpha_cens,
      scaling = scaling,
      BG = BG,
      nsim = final_n + 2L*Nk,
      posterior_mode = "dataset",
      posterior_draw_index = draw_j
    )

    n0 <- sum(baseline_dec$treat == 0)
    n1 <- sum(baseline_dec$treat == 1)
    if (min(n0,n1) <= Nk) {
      stop("Too few subjects in one arm for fixed-n decoupling at N=", Nk,
           ". n0=",n0,", n1=",n1)
    }

    d_dec <- decoupling_covariate_surv1(
      N = Nk,
      p = 0.5,
      x = x_decoupling,
      beta_weib = weib$beta_event,
      beta_weib_cens = weib$beta_cens,
      sim_data = baseline_dec
    )
    if (nrow(d_dec) != final_n) stop("Decoupling fixed-n construction failed.")

    one <- GBSG2_diag_one(d_dec, source_data = datt, mahal_index = c(1,9,10))
    rr <- as.data.frame(as.list(one))
    rr$method <- "Decoupling Covariate"
    rr$N <- Nk
    rr$replicate <- j
    rr$final_n_target <- final_n
    rep_rows[[zz]] <- rr
    zz <- zz + 1L
  }
}

replicates <- do.call(rbind, rep_rows)

summary_table <- do.call(rbind, lapply(
  split(replicates, interaction(replicates$method, replicates$N, drop = TRUE)),
  function(m) {
    s <- GBSG2_summarise_replicates(
      m,
      extra = list(
        method = m$method[1],
        N = m$N[1],
        final_n_target = final_n
      )
    )
    s[, c("method","N","final_n_target",
          setdiff(names(s), c("method","N","final_n_target")))]
  }
))
rownames(summary_table) <- NULL

write.csv(summary_table,
          "results/GBSG2_fixed_final_n_summary.csv",
          row.names = FALSE)
write.csv(replicates,
          "results/GBSG2_fixed_final_n_replicates.csv",
          row.names = FALSE)

cat("\nGBSG2 fixed-final-n sensitivity:\n")
print(summary_table, row.names = FALSE)

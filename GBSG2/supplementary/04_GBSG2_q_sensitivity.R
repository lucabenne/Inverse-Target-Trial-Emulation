###############################################################################
# 04_GBSG2_q_sensitivity.R
#
# Sensitivity of Additional Sampling to the source seed-pair quantile q.
# The final v6 main setting q=0.80 is included.
# N is fixed at the maximum main-analysis level (N=100).
###############################################################################

if (!exists("weib", inherits = FALSE) || !exists("BG", inherits = FALSE)) {
  
}
source("99_GBSG2_diagnostic_helpers.R")

q_levels <- c(0.50, 0.65, 0.80, 0.95)
N_add <- 100L
n_rep <- 1000L

index_additional <- c(1,2,3,9,10)
Ma_source <- pairwise_mahalanobis(datt, index_additional)
Ma_sorted <- sort(Ma_source)

seed_distance <- vapply(q_levels, function(q) {
  pos <- max(1L, min(length(Ma_sorted), round(length(Ma_sorted) * q)))
  Ma_sorted[pos]
}, numeric(1))

set.seed(140501L)
BG_draws <- sample.int(n_BG_draws, n_rep, replace = TRUE)

rep_rows <- list()
zz <- 1L

for (j in seq_len(n_rep)) {
  # Same Bernoulli addition sequence for every q within the replicate.
  set.seed(141000L + j)
  bern_j <- rbinom(N_add, 1, 0.5)

  for (k in seq_along(q_levels)) {
    qk <- q_levels[k]

    # Reset the stochastic local generation for each q to a deterministic,
    # q-specific stream while leaving the underlying mechanism unchanged.
    set.seed(142000L + j*100L + k)

    d <- newdatamahal_distr_surv2(
      data = datt,
      bern = bern_j,
      n = N_add,
      q = qk,
      beta_weib = weib$beta_event,
      beta_weib_cens = weib$beta_cens,
      alpha_event = weib$alpha_event,
      alpha_cens = weib$alpha_cens,
      scaling = scaling,
      Ma = Ma_source,
      BG = BG,
      data_prec = datt,
      BG_draw_index = BG_draws[j]
    )

    one <- GBSG2_diag_one(d, source_data = datt, mahal_index = index_additional)
    rr <- as.data.frame(as.list(one))
    rr$q <- qk
    rr$seed_pair_mahalanobis <- seed_distance[k]
    rr$N <- N_add
    rr$replicate <- j
    rep_rows[[zz]] <- rr
    zz <- zz + 1L
  }
}

replicates <- do.call(rbind, rep_rows)

summary_table <- do.call(rbind, lapply(
  split(replicates, replicates$q),
  function(m) {
    s <- GBSG2_summarise_replicates(
      m,
      extra = list(
        q = m$q[1],
        seed_pair_mahalanobis = m$seed_pair_mahalanobis[1],
        N = N_add
      )
    )
    s[, c("q","seed_pair_mahalanobis","N",
          setdiff(names(s), c("q","seed_pair_mahalanobis","N")))]
  }
))
rownames(summary_table) <- NULL

write.csv(summary_table,
          "results/GBSG2_q_sensitivity_summary.csv",
          row.names = FALSE)
write.csv(replicates,
          "results/GBSG2_q_sensitivity_replicates.csv",
          row.names = FALSE)

cat("\nGBSG2 q-sensitivity summary:\n")
print(summary_table, row.names = FALSE)

###############################################################################
# 02_GBSG2_severity_overlap_diagnostics.R
#
# Reviewer-requested diagnostics beyond Mahalanobis:
# SMD, PS overlap, ESS / extreme weights, and pooled source-vs-final shift.
#
# Uses the frozen v6 GBSG2 mechanism settings:
# N = 20,40,60,80,100; q=0.80; 1,000 Monte Carlo replicates.
###############################################################################

if (!exists("weib", inherits = FALSE) || !exists("BG", inherits = FALSE)) {
  
}
source("99_GBSG2_diagnostic_helpers.R")

N_levels <- c(20L,40L,60L,80L,100L)
n_rep <- 1000L
q_additional <- 0.80
x_decoupling <- "pnodes"
nsim <- 1000L

index_additional <- c(1,2,3,9,10)
Ma_source <- pairwise_mahalanobis(datt, index_additional)

set.seed(100501L)
BG_draw_additional <- sample.int(n_BG_draws, n_rep, replace = TRUE)
BG_draw_decoupling <- sample.int(n_BG_draws, n_rep, replace = TRUE)

rep_rows <- list()
smd_rows <- list()
zz <- 1L

# Additional Sampling: cumulative N levels exactly as in the final v6 main code.
for (j in seq_len(n_rep)) {
  set.seed(110000L + j)
  bern_j <- rbinom(max(N_levels), 1, 0.5)
  data_add <- datt
  previous_N <- 0L

  for (k in seq_along(N_levels)) {
    n_add <- N_levels[k] - previous_N
    idx_bern <- seq.int(previous_N + 1L, N_levels[k])

    data_add <- newdatamahal_distr_surv2(
      data = datt,
      bern = bern_j[idx_bern],
      n = n_add,
      q = q_additional,
      beta_weib = weib$beta_event,
      beta_weib_cens = weib$beta_cens,
      alpha_event = weib$alpha_event,
      alpha_cens = weib$alpha_cens,
      scaling = scaling,
      Ma = Ma_source,
      BG = BG,
      data_prec = data_add,
      BG_draw_index = BG_draw_additional[j]
    )

    one <- GBSG2_diag_one(data_add, source_data = datt, mahal_index = index_additional)
    rr <- as.data.frame(as.list(one))
    rr$method <- "Additional Sampling"
    rr$N <- N_levels[k]
    rr$replicate <- j
    rep_rows[[zz]] <- rr

    s <- attr(one, "smd")
    smd_rows[[zz]] <- data.frame(
      method = "Additional Sampling",
      N = N_levels[k],
      replicate = j,
      covariate = names(s),
      smd = as.numeric(s),
      abs_smd = abs(as.numeric(s)),
      row.names = NULL
    )
    zz <- zz + 1L
    previous_N <- N_levels[k]
  }
}

# Decoupling Covariate: one fresh synthetic RCT per replicate, same RCT across N.
for (j in seq_len(n_rep)) {
  set.seed(120000L + j)
  sim_data <- generate_trial(
    data = datt,
    beta_weib = weib$beta_event,
    beta_weib_cens = weib$beta_cens,
    alpha_event = weib$alpha_event,
    alpha_cens = weib$alpha_cens,
    scaling = scaling,
    BG = BG,
    nsim = nsim,
    posterior_mode = "dataset",
    posterior_draw_index = BG_draw_decoupling[j]
  )

  for (k in seq_along(N_levels)) {
    d <- decoupling_covariate_surv1(
      N = N_levels[k],
      p = 0.5,
      x = x_decoupling,
      beta_weib = weib$beta_event,
      beta_weib_cens = weib$beta_cens,
      sim_data = sim_data
    )

    one <- GBSG2_diag_one(d, source_data = datt, mahal_index = c(1,9,10))
    rr <- as.data.frame(as.list(one))
    rr$method <- "Decoupling Covariate"
    rr$N <- N_levels[k]
    rr$replicate <- j
    rep_rows[[zz]] <- rr

    s <- attr(one, "smd")
    smd_rows[[zz]] <- data.frame(
      method = "Decoupling Covariate",
      N = N_levels[k],
      replicate = j,
      covariate = names(s),
      smd = as.numeric(s),
      abs_smd = abs(as.numeric(s)),
      row.names = NULL
    )
    zz <- zz + 1L
  }
}

replicates <- do.call(rbind, rep_rows)
smd_by_covariate <- do.call(rbind, smd_rows)

summary_table <- do.call(rbind, lapply(
  split(replicates, interaction(replicates$method, replicates$N, drop = TRUE)),
  function(m) {
    s <- GBSG2_summarise_replicates(
      m,
      extra = list(method = m$method[1], N = m$N[1])
    )
    s[, c("method","N", setdiff(names(s), c("method","N")))]
  }
))
rownames(summary_table) <- NULL

write.csv(summary_table,
          "results/GBSG2_severity_overlap_summary.csv",
          row.names = FALSE)
write.csv(replicates,
          "results/GBSG2_severity_overlap_replicates.csv",
          row.names = FALSE)
write.csv(smd_by_covariate,
          "results/GBSG2_severity_SMD_by_covariate.csv",
          row.names = FALSE)

cat("\nGBSG2 severity / overlap summary:\n")
print(summary_table, row.names = FALSE)

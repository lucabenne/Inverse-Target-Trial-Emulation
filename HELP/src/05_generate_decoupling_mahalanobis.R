decoupling_mahalanobis2 = function(data,N,p,x1,resampling,n_res,BG){

  library(readr)
  library(MatchIt)

  if (!exists("generate_HELP_synthetic_trial")) {
    source("03a_generate_HELP_synthetic_trial.R")
  }

  # Each Monte Carlo replicate starts from a fresh randomized synthetic HELP RCT.
  # The Bayesian DGP is fitted once upstream; baseline-covariate parameters use
  # one joint posterior draw per replicate, while outcome parameters are fixed
  # at posterior means so the treatment-effect target stays constant.

  if (length(N) == 1L && N == 0) {
    return(generate_HELP_synthetic_trial(
      data = data,
      BG = BG,
      nsim = 1000L,
      treat_prob = p,
      posterior_mode = "dataset",
      outcome_posterior_mode = "mean"
    ))
  }

  make_one = function(N_target) {

    data_rep = generate_HELP_synthetic_trial(
      data = data,
      BG = BG,
      nsim = 1000L,
      treat_prob = p,
      posterior_mode = "dataset",
      outcome_posterior_mode = "mean"
    )

    target_rct = attr(data_rep, "rct_target")

    if (N_target == 0) return(data_rep)

    if (resampling == 1) {

      data_cov = data_rep[,-c(1,3)]
      ind_nottreat = which(data_rep$treat == 0)
      ind_treat = which(data_rep$treat == 1)

      treat_cov = data_cov[ind_treat,,drop=FALSE]
      nottreat_cov = data_cov[ind_nottreat,,drop=FALSE]

      ncol_cov = ncol(data_cov)
      length_nottreat = length(ind_nottreat)
      length_treat = length(ind_treat)

      Ma = matrix(0,length_nottreat,length_treat)
      V = cov(data_cov)

      for (j in seq_len(length_treat)) {
        for (i in seq_len(length_nottreat)) {
          d = as.numeric(treat_cov[j,] - nottreat_cov[i,])
          e = solve(V,d)
          Ma[i,j] = sum(d*e)
        }
      }

      x_treat = data.frame(
        index = seq_along(ind_treat),
        x = data_rep[[x1]][ind_treat]
      )
      x_treat = x_treat[order(x_treat$x),]

      sampled_individual = numeric()
      index_to_eliminate = numeric()
      i = 1L
      j = 1L

      while (length(index_to_eliminate) < N_target) {
        sample_1_n = runif(1,1,n_res)
        sorted_distance = sort(Ma[,x_treat$index[i]])
        rank_id = min(round(sample_1_n),length(sorted_distance))
        sampled_individual[i] = which(Ma[,x_treat$index[i]] == sorted_distance[rank_id])[1]
        index_to_eliminate[i] = ind_nottreat[sampled_individual[i]]

        while (length(which(table(index_to_eliminate)>1)) > 0) {
          rank_id = min(round(sample_1_n)+j,length(sorted_distance))
          sampled_individual[i] = which(Ma[,x_treat$index[i]] == sorted_distance[rank_id])[1]
          index_to_eliminate[i] = ind_nottreat[sampled_individual[i]]
          j = j+1L
        }

        i = i+1L
        j = 1L
      }

      out = data_rep[-index_to_eliminate,,drop=FALSE]
      rownames(out) = NULL
      attr(out, "rct_target") = target_rct
      return(out)
    }

    # Original no-resampling Decoupling Mahalanobis mechanism.
    data_rep = data_rep[order(data_rep[[x1]]),]
    rownames(data_rep) = seq_len(nrow(data_rep))

    formula <- as.formula(
      "treat ~ female + daysdrink + daysanysub + pcs + mcs + cesd + sexrisk + drugrisk"
    )

    m.out = matchit(
      formula = formula,
      data = data_rep,
      distance = "mahalanobis"
    )

    x = as.numeric(m.out$match.matrix[,1])
    y = data.frame(
      treat = m.out$treat,
      index = seq_len(nrow(data_rep))
    )
    y = y[y$treat == 1,,drop=FALSE]
    x = data.frame(
      treat = y$treat,
      V2 = y$index,
      x = x
    )[,c("V2","x")]

    ind_treat = which(data_rep$treat == 1)
    x_treat = data.frame(
      index = ind_treat,
      x = data_rep[[x1]][ind_treat]
    )

    removed_individuals = numeric()
    i = 1L
    j = 1L

    while (i <= N_target) {

      if ((i+j) > nrow(x_treat)) {
        stop("Not enough matched treated subjects to remove the requested number of controls.")
      }

      index_to_eliminate = x[x$V2 == x_treat$index[i],"x"]

      while (length(index_to_eliminate) == 0L || is.na(index_to_eliminate)) {
        index_to_eliminate = x[x$V2 == x_treat$index[i+j],"x"]
        j = j+1L
        if ((i+j) > nrow(x_treat)+1L &&
            (length(index_to_eliminate) == 0L || is.na(index_to_eliminate))) {
          stop("Unable to identify enough matched controls for decoupling.")
        }
      }

      removed_individuals[i] = index_to_eliminate
      i = i+1L
      j = 1L
    }

    out = data_rep[-removed_individuals,,drop=FALSE]
    rownames(out) = NULL
    attr(out, "rct_target") = target_rct
    out
  }

  if (length(N) == 1L) return(make_one(N))

  list_data = vector("list",length(N))
  for (k in seq_along(N)) {
    if (k == 1L || k %% 25L == 0L || k == length(N)) {
      cat("Decoupling Mahalanobis replicate",k,"of",length(N),"\n")
    }
    list_data[[k]] = make_one(N[k])
  }
  list_data
}

pstdm <- function(sim_data){
  # PTDM mechanism unchanged; indexing written explicitly for robustness.
  sim_data1 <- sim_data
  index_treated <- which(sim_data1$new_treat==1)
  index_nottreated <- which(sim_data1$new_treat==0)

  if(!length(index_treated)) stop("PTDM requires at least one treated individual.")
  if(length(index_nottreated)){
    sim_data1$immortal_time[index_nottreated] <- sample(
      sim_data1$immortal_time[index_treated],length(index_nottreated),replace=TRUE
    )
  }

  index_to_exclude <- which(sim_data1$immortal_time>sim_data1$time)
  if(length(index_to_exclude)) sim_data1 <- sim_data1[-index_to_exclude,,drop=FALSE]

  sim_data1$time <- sim_data1$time-sim_data1$immortal_time
  sim_data1$treat <- sim_data1$new_treat
  sim_data1
}

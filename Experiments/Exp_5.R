#### Experiments, Nov 15th
#### Examining community recovery as the number of nodes is varied

.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)

source(here("Experiments/", "utils.R"))

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

node_vec <- c(50, 100, 200, 500, 1000)

no_sims <- 20
curr_n <- node_vec[sim_id]


results <- list()

for(sim in 1:no_sims) {
  cat("Sim:", sim, "\n")
  n <- curr_n
  Time <- 100
  intens1 <- c(2)
  intens2 <- c(1)
  true_B <- matrix(c(intens1, 0.05, 0.05, intens2), 
                   nrow = 2, ncol = 2, byrow = T)
  
  intens <- matrix(c(intens1, 0.05, 0.05, intens2), 4, 1)
  sim1 <- generateDynppsbmConst(intens = intens,
                                Time = Time,
                                n = n, 
                                prop.groups = c(0.5, 0.5))
  proc_sim <- format_sims(sim_data = sim1, n = n)
  
  K <- 2
  m <- n
  B <- matrix(runif(K * K), K, K)
  dT <- curr_dT #dT_vec[sim] # this should be current dT?
  cat("dT:", curr_dT, "\n")
  inter_T <- 1
  
  results_online <- estimate_Poisson(full_data = proc_sim$events,
                                     A = proc_sim$edge,
                                     m,
                                     K,
                                     Time,
                                     dT = curr_dT,
                                     B,
                                     inter_T,
                                     is_elbo = TRUE)
  
  # compute rand index
  z_true <- apply(sim1$z, 2, function(x) which(x == 1))
  z_est <- apply(results_online$tau, 1, which.max)
  clust_est <- aricode::ARI(z_true, z_est)
  
  sim_pars <- list(
    B = B, 
    est_elbo = results_online$AveELBO,
    clust = clust_est,
    tau = results_online$tau,
    curr_n = n
  )
  results[[sim]] <- sim_pars
}

# all_results <- list(results, init)

saveRDS(results, file = here("Experiments",
                             "exp_results",
                             paste0("exp5_", curr_n, ".RDS")))

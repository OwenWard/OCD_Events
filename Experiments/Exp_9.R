#### Experiments, Dec 6th
#### To explore the online community recovery for simulation models
.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)

source(here("Experiments/", "utils.R"))

nsims <- 20

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid
Times <- c(50, 100, 200, 500)
Time <- Times[sim_id]

results <- list()

for(sim in 1:nsims) {
  cat("Sim:", sim, "\n")
  # Time <- 100
  n <- 100
  intens1 <- c(3)
  intens2 <- c(1)
  intens <- matrix(c(intens1, 0.05, 0.05, intens2), 4, 1)
  true_B <- matrix(c(intens1, 0.05, 0.05, intens2), 
                   nrow = 2, ncol = 2, byrow = T)
  # this is essentially the K*K matrix stretched out as a single col
  system.time(sim1 <- generateDynppsbmConst(intens = intens,
                                            Time = Time,
                                            n = n, 
                                            prop.groups = c(0.5, 0.5)))
  
  proc_sim <- format_sims(sim_data = sim1, n = n)
  K <- 2
  m <- n
  B <- matrix(runif(K * K), K, K)
  dT <- 1
  inter_T <- 1
  # capture output to not print out
  results_online <- estimate_Poisson(full_data = proc_sim$events,
                                     A = proc_sim$edge,
                                     m,
                                     K,
                                     Time,
                                     dT,
                                     B,
                                     inter_T,
                                     is_elbo = TRUE)
  
  
  B_ests <- results_online$inter_B
  z_true <- apply(sim1$z, 2, which.max)
  z_est <- apply(results_online$tau, 1, which.max)
  clust_est <- aricode::ARI(z_true, z_est)
  
  sim_pars <- list(
    B = B,
    est_elbo = results_online$AveELBO,
    clust = clust_est,
    online_clust = results_online$early_tau,
    z_true = z_true
  )
  results[[sim]] <- sim_pars
}

saveRDS(results, file = here("Experiments",
                             "exp_results",
                             paste0("exp9_", Time, ".RDS")))


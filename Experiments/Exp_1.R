#### Experiments, Oct 26th
#### To investigate the relationship between the 
#### the initial B matrix and the clustering result
#### and confirm there is a scenario where performance is good
#### for the majority of initial values, and see how it improves with time
.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)

source(here("Experiments/", "utils.R"))

nsims <- 100
all_times <- c(50, 100, 250, 500, 1000, 2000)

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

results <- list()
Time <- all_times[sim_id]

for(sim in 1:nsims) {
  cat("Sim:", sim, "\n")
  n <- 100
  # intens1 <- c(2.2)
  # intens2 <- c(1.75)
  # true_B <- matrix(c(intens1, 0.35, 0.5, intens2), 
  #                  nrow = 2, ncol = 2, byrow = T)
  # 
  # intens <- matrix(c(intens1, 0, 0.1, intens2), 4, 1)
  intens1 <- c(2)
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
  # compute rand index
  z_true <- apply(sim1$z, 2, which.max)
  z_est <- apply(results_online$tau, 1, which.max)
  clust_est <- aricode::ARI(z_true, z_est)
  
  sim_pars <- list(
    B = B, 
    est_B = results_online$B
    est_elbo = results_online$AveELBO,
    clust = clust_est,
    time = Time
  )
  results[[sim]] <- sim_pars
}

saveRDS(results, file = here("Experiments",
                             "exp_results",
                             paste0("exp1_store_B_time_", Time, ".RDS")))

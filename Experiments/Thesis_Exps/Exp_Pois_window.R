#### Experiment Pois Window, May 10th
#### to investigate the relationship between 
#### the window size and the clustering result,
#### for a fixed data set and initial B matrix

.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)

source(here("Experiments/", "utils.R"))
source(here("functions/init_fcn.R"))

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

dT_vec <- seq(from = 0.25, to = 5, by = 0.25)
## 50 scenarios here now

sparsity <- 0.05
no_sims <- 50
curr_dT <- dT_vec[sim_id]


results <- tibble()

for(sim in 1:no_sims) {
  cat("Sim:", sim, "\n")
  m <- 200
  Time <- 200
  m0_curr <- m/2
  curr_n0 <- 20
  inter_T <- 1
  K <- 3
  ## baseline rate of the process
  true_Mu <- matrix(c(2, 0.05, 0.15, 0.05, 1.25, 0.05, 0.15, 0.05, 0.75), 
                    nrow = K, ncol = K, byrow = T)
  
  true_B <- matrix(0, nrow = K, ncol = K, byrow = TRUE)
  # Pi <- c(0.2, 0.3, 0.3, 0.2)
  Pi <- c(1/K, 1/K, 1/K)
  
  Z <- sample(0:(K-1), size = m, prob = Pi, replace = TRUE)
  # then generate A
  A <- list()
  for(i in 1:m){
    # could sample these here with SBM structure...
    num_edge = m * sparsity
    edge <- sample(m, num_edge) - 1
    edge <- sort(edge[edge!= i-1]) ## remove self edges
    A[[i]] <- edge
  }
  alltimes <- sampleBlockHak(Time, A, Z, Mu = true_Mu, B = true_B, lam = 1)
  
  result <- sparse_poisson(alltimes, K, n0 = curr_n0, m, m0 = m0_curr)

  
  Mu_est <- result$est_B
  ## need to pass the estimated clustering also
  init_tau <- matrix(0, nrow = m, ncol = K)
  for(i in seq_along(result$est_clust)){
    init_tau[i, result$est_clust[i]] <- 1
  }
  
  results_online_init <- estimate_Poisson_init(full_data = 
                                                 result$rest_events,
                                               A,
                                               m,
                                               K,
                                               Time,
                                               dT = curr_dT,
                                               B = Mu_est,
                                               inter_T,
                                               init_tau,
                                               start = result$cut_off,
                                               is_elbo = FALSE)
  z_est <- apply(results_online_init$tau, 1, which.max)
  clust_est_init <- aricode::ARI(Z, z_est)
  

  
  sim_pars <- tibble(
    # B = Mu_est, 
    sim = sim,
    clust = clust_est_init,
    # tau = results_online_init$tau,
    dT = curr_dT
  )
  results <- results %>%  bind_rows(sim_pars)
}

# all_results <- list(results, init)

saveRDS(results, file = here("Experiments",
                             "thesis_output",
                             paste0("exp_pois_wind_may_23_", curr_dT, ".RDS")))
# to save it after each run also


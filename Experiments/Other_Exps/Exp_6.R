#### Experiment 6, Nov 15th
#### Examining community recovery as the number of groups is varied

.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)

source(here("Experiments/", "utils.R"))

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

k_vec <- c(3:6)

no_sims <- 50
curr_k <- k_vec[sim_id]

### need to form the intensity matrix here in a smart way
form_intensity <- function(k, within_rates, between_rate) {
  stopifnot(length(within_rates) == k)
  B <- matrix(rep(between_rate, k*k), nrow = k, ncol = k)
  diag(B) <- within_rates
  return(as.vector(t(B)))
}

all_rates <- c(1:6)
cur_rates <- all_rates[1:curr_k]


results <- list()

for(sim in 1:no_sims) {
  cat("Sim:", sim, "\n")
  n <- 100
  Time <- 100
  intens1 <- c(2)
  intens2 <- c(1)
  true_B <- matrix(c(intens1, 0.05, 0.05, intens2), 
                   nrow = 2, ncol = 2, byrow = T)
  
  # intens <- matrix(c(intens1, 0.05, 0.05, intens2), 4, 1)
  intens <- matrix(form_intensity(curr_k, cur_rates, 0.05), curr_k*curr_k, 1)
  sim1 <- generateDynppsbmConst(intens = intens,
                                Time = Time,
                                n = n, 
                                prop.groups = rep(1/curr_k, curr_k))
  proc_sim <- format_sims(sim_data = sim1, n = n)
  
  K <- curr_k
  m <- n
  B <- matrix(runif(K * K), K, K)
  dT <- 1 #dT_vec[sim] # this should be current dT?
  # cat("dT:", curr_dT, "\n")
  inter_T <- 1
  
  results_online <- estimate_Poisson(full_data = proc_sim$events,
                                     A = proc_sim$edge,
                                     m,
                                     K,
                                     Time,
                                     dT = dT,
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
    K = curr_k
  )
  results[[sim]] <- sim_pars
}

# all_results <- list(results, init)

saveRDS(results, file = here("Experiments",
                             "exp_results",
                             paste0("exp6_", curr_k, ".RDS")))

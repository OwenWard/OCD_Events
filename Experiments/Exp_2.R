#### Experiment 2, Oct 26th
#### to investigate the relationship between 
#### the window size and the clustering result,
#### for a fixed data set and initial B matrix

.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)

source(here("Experiments/", "utils.R"))

# jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# jobid <- as.numeric(jobid)
# sim_id <- jobid

dT_vec <- seq(from = 0.1, to = 2, by = 0.1)

no_sims <- length(dT_vec)

### then load in a dataset and an initial B? once which works 
### well for a known dT
exp <- readRDS(here("Experiments/", "sim_pars.RDS"))

B <- exp$init_B
sim1 <- exp$sim
n <- 100

proc_sim <- format_sims(sim_data = sim1, n = n)

K <- 2
m <- n
Pi <- rep(1/K, K)
tau <- matrix(1, nrow = m, ncol = K)
tau <- tau/rowSums(tau)
S <- matrix(1/K, nrow = m, ncol = K)


inter_T <- 1
results <- list()

for(sim in 1:no_sims) {
  dT <- dT_vec[sim]
  
  results_online <- estimate_Poisson(full_data = proc_sim$events,
                                     A = proc_sim$edge,
                                     m,
                                     K,
                                     Time,
                                     dT,
                                     B,
                                     tau,
                                     Pi,
                                     S,
                                     inter_T,
                                     is_elbo = TRUE)
  
  # compute rand index
  z_true <- apply(sim1$z, 2, function(x) which(x == 1))
  z_est <- apply(results_online$tau, 1, which.max)
  clust_est <- aricode::ARI(z_true, z_est)
  
  sim_pars <- list(
    dT = dT, 
    est_eblo = results_online$AveELBO,
    clust = clust_est
  )
  results[[sim]] <- sim_pars
}


saveRDS(results, file = here("Experiments", "exp_dt.RDS"))
#### Exp verifying online community recovery, Apr 26th
#### Investigate whether the initialization function 
#### leads to an improved clustering performance
.libPaths("/moto/stats/users/ogw2103/rpackages")
library(here)

source(here("Experiments/", "utils.R"))
source(here("functions/init_fcn.R"))


### simulate some data, see how the performance changes with the
### init function


Time <- 200
no_sims <- 50
dT <- 1
inter_T <- 20
K <- 2
m_vec <- c(200, 500, 1000)

sparsity <- 0.1 # prop of edges which can have events

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

model <- "Hawkes"

results <- list()
m <- m_vec[sim_id]


# for(exp_num in seq_along(m_vec)) {
#   dT <- 1
curr_dt_sims <- tibble()
cat("Current K:", K, "\n")
cat("Current m:", m, "\n")
cat(model, "\n")
for(sim in 1:no_sims){
  cat("Sim:", sim, "======\n")
  ## baseline rate of the process
  true_Mu <- matrix(0.25, 
                    nrow = K, ncol = K, byrow = T)
  diag(true_Mu) <- 1:K
  true_B <- matrix(0, nrow = K, ncol = K, byrow = TRUE)
  diag(true_B) <- 0.5
  # Pi <- c(0.2, 0.3, 0.3, 0.2)
  Pi <- c(0.2, 0.3)
  
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
  
  ## random initialization for now
  Mu <- matrix(runif(K * K), K, K)
  B <- matrix(runif(K * K), K, K)
  tau <- matrix(1/K, nrow = m, ncol = K)
  results_online <- online_estimator_eff_revised(alltimes, 
                                                 A,
                                                 m,
                                                 K,
                                                 Time,
                                                 dT = dT,
                                                 lam = 1,
                                                 B, 
                                                 Mu,
                                                 tau,
                                                 inter_T, 
                                                 is_elbo = FALSE)
  
  ### compute the online community assignments here
  # init_ari <- aricode::ARI(result$est_clust, Z)
  inter_tau <- results_online$early_tau[, , 1:10]
  inter_z <- apply(inter_tau, c(1, 3), which.max)
  inter_ari <- apply(inter_z, 2, function(x) aricode::ARI(x, Z))
  z_est <- apply(results_online$tau, 1, which.max)
  ari_final <- aricode::ARI(Z, z_est)
  ###
  ### then save dT, clust_est, m, model
  curr_sim <- tibble(ARI = inter_ari,
                     K = K, 
                     nodes = m,
                     model = model,
                     time = seq(from = inter_T, to = Time, by = inter_T),
                     sparsity = sparsity,
                     sim = sim)
  curr_dt_sims <- curr_dt_sims %>% 
    bind_rows(curr_sim)

}



results <- curr_dt_sims

### then save these somewhere
saveRDS(results, file = here("Experiments",
                             "thesis_output",
                             paste0("exp_hawkes_new_online_rho_",
                                    100*sparsity, sim_id, ".RDS")))

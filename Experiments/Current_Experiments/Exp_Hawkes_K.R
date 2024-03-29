#### Exp varying the number of Nodes, April 25th 2022
#### Investigate whether the initialization function 
#### leads to an improved clustering performance
.libPaths("/moto/stats/users/ogw2103/rpackages")
library(here)

source(here("functions/utils.R"))
source(here("functions/init_fcn.R"))
source(here("functions/init_fcn_Hawkes.R"))

### simulate some data, see how the performance changes with the
### init function


Time <- 200
no_sims <- 50
dT <- 1
inter_T <- 1
K_vec <- c(3, 4, 5, 6)
# m_vec <- c(rep(100, 3), rep(200, 3), rep(500, 3), rep(1000, 3),
#            rep(5000, 3))
m <- 200
sparsity <- 0.15 # prop of edges which can have events

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

model <- "Hawkes"
# if(sim_id <= 3){
#   model = "Poisson"
# }

results <- list()

K <- K_vec[sim_id]


m0 <- m/4
n0 <- 20

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
  diag(true_Mu) <- c(1:K)
  true_B <- matrix(0.0, nrow = K, ncol = K, byrow = TRUE)
  diag(true_B) <- 0.5
  # Pi <- c(0.2, 0.3, 0.3, 0.2)
  Pi <- rep(1/K, K)
  
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
  print("Simulated Data")
  result <- sparse_Hawkes(alltimes,
                          K,
                          n0 = n0,
                          m, m0)
  Mu_init <- result$est_B[, , 1]
  B_init <- result$est_B[, , 2]
  ## need to pass the estimated clustering also
  init_tau <- matrix(0, nrow = m, ncol = K)
  for(i in seq_along(result$est_clust)){
    init_tau[i, result$est_clust[i]] <- 1
  }
  S <- matrix(1/K, nrow = m, ncol = K)
  results_online <- online_estimator_eff_revised(alltimes, 
                                                 A,
                                                 m,
                                                 K,
                                                 Time,
                                                 dT = dT,
                                                 lam = 1,
                                                 B_init, 
                                                 Mu_init,
                                                 init_tau,
                                                 S,
                                                 inter_T, 
                                                 is_elbo = FALSE)
  
  
  stan_est <- apply(results_online$tau, 1, which.max)
  (clust_est_norm <- aricode::ARI(stan_est, Z))
  curr_sim_rand <- tibble(ARI = clust_est_norm,
                          K = K, 
                          nodes = m,
                          model = model,
                          init = "Init",
                          n0 = n0,
                          m0 = m0,
                          sparsity = sparsity,
                          sim = sim_id)
  curr_dt_sims <- curr_dt_sims %>% 
    bind_rows(curr_sim_rand)
  
}

results <- curr_dt_sims

### then save these somewhere
saveRDS(results, file = here("Experiments",
                             "exp_results", "Sept_23",
                             paste0("exp_hawkes_k_",
                                    K, "_", sim_id, ".RDS")))


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
K <- 2
# m_vec <- c(rep(100, 3), rep(200, 3), rep(500, 3), rep(1000, 3),
#            rep(5000, 3))
# m_vec <- c(rep(100, 3), rep(200, 3), rep(500, 3), rep(1000, 3), rep(5000, 3))
m_vec <- c(100, 200, 500, 1000)
m <- 1000

sparsity <- 0.15 # prop of edges which can have events

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

model <- "Hawkes"

results <- list()
m <- m_vec[sim_id]
# m0_vec <- c( 100*c(1/10, 1/4, 1/2),
#              200*c(1/10, 1/4, 1/2),
#              500*c(1/10, 1/4, 1/2),
#              1000*c(1/10, 1/4, 1/2))
# m0_vec <- c(10, 20, 50, 100)
# m0_curr <- m0_vec[sim_id]

if(m == 1000){
  no_sims <- 5
}


m0 <- m/4
n0 <- 40
inter_T <- 1

# for(exp_num in seq_along(m_vec)) {
#   dT <- 1
init_results <- tibble()

curr_dt_sims <- tibble()
cat("Current K:", K, "\n")
cat("Current m:", m, "\n")
cat(model, "\n")
for(sim in 1:no_sims){
  cat("Sim:", sim, "======\n")
  ## baseline rate of the process
  true_Mu <- matrix(c(1, 0.25, 0.25, 2), 
                    nrow = K, ncol = K, byrow = T)
  true_B <- matrix(c(0.5, 0, 0, .5), nrow = K, ncol = K, byrow = TRUE)
  # diag(true_B) <- 0.5
  if(model == "Poisson") {
    true_B <- matrix(0, K, K)
  }
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
  print("Simulated Data")
  
  result <- sparse_Hawkes(alltimes,
                          K,
                          n0 = n0,
                          m, m0)
  Mu_est <- result$est_B
  ## need to pass the estimated clustering also
  init_tau <- matrix(0, nrow = m, ncol = K)
  for(i in seq_along(result$est_clust)){
    init_tau[i, result$est_clust[i]] <- 1
  }
  ### check the initial ARI
  cat("Init Scheme \n")
  print(aricode::ARI(result$est_clust, Z))
  init_results <- init_results %>% bind_rows(tibble(init = "Scheme",
                                    ARI = aricode::ARI(result$est_clust, Z)))
  cat("-------\n")
  init_tau <- matrix(0, nrow = m, ncol = K)
  for(i in seq_along(result$est_clust)){
    init_tau[i, result$est_clust[i]] <- 1
  }
  Mu_init <- result$est_B[, , 1]
  B_init <- result$est_B[, , 2]
  ## random initialization for now
  Mu <- matrix(runif(K * K), K, K)
  B <- matrix(runif(K * K), K, K)
  tau <- matrix(1/K, nrow = m, ncol = K)
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
  cat("Random \n")
  (clust_est_norm <- aricode::ARI(stan_est, Z))
  print(clust_est_norm)
  # init_results <- init_results %>% 
  #   bind_rows(tibble(init = "Random", ARI = clust_est_norm))
  cat("-----\n")
  # print("Random Worked")
  curr_sim_rand <- tibble(ARI = clust_est_norm,
                          K = K, 
                          nodes = m,
                          model = model,
                          init = "Init",
                          n0 = NA,
                          m0 = NA,
                          sparsity = sparsity,
                          sim = sim_id)
  curr_dt_sims <- curr_dt_sims %>% 
    bind_rows(curr_sim_rand)
  
}

# init_results_curr <- init_results
# init_results_prev <- init_results

results <- curr_dt_sims

### then save these somewhere
if(m < 1000){
  saveRDS(results, file = here("Experiments",
                               "exp_results", "Sept_23",
                               paste0("exp_hawkes_nodes_",
                                      m, "_", sim_id, ".RDS")))
}
if(m == 1000){
  saveRDS(results, file = here("Experiments",
                               "exp_results", "Sept_23",
                               paste0("exp_hawkes_nodes_1k_",
                                      m, "_", sim_id, ".RDS")))
}

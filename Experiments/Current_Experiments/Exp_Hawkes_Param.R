#### Exp verifying online community recovery, Apr 26th
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
inter_T <- 20
K <- 2
Time_vec <- c(50, 100, 200, 500)
m <- 200

m0 <- m/4
n0 <- 20

sparsity <- 0.15 # prop of edges which can have events

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

model <- "Hawkes"

results <- list()
##m <- m_vec[sim_id]
Time <- Time_vec[sim_id]

# for(exp_num in seq_along(m_vec)) {
#   dT <- 1
curr_dt_sims <- tibble()
cat("Current K:", K, "\n")
cat("Current m:", m, "\n")
cat(model, "\n")
for(sim in 1:no_sims){
  cat("Sim:", sim, "======\n")
  ## baseline rate of the process
  true_Mu <- matrix(c(2, 0.05, 0.1, 1.5), 
                    nrow = K, ncol = K, byrow = T)
  
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
  
  z_est <- apply(results_online$tau, 1, which.max)
  ari_final <- aricode::ARI(Z, z_est)
  
  sim_pars <- list(
    z_true = Z,
    z_est = z_est,
    ## and these
    B_ests = results_online$inter_B,
    Mu_ests = results_online$inter_mu,
    clust = ari_final,
    lambda = results_online$lambda_vec,
    # regret = regret,
    # card_A = card_A,
    # batch_ave_loss = batch_average,
    # online_loss = tidy_loss,
    # pred_llh = ave_pred_ll,
    Time = Time
  )
  results[[sim]] <- sim_pars
  
}

# results 

### then save these somewhere
saveRDS(results, file = here("Experiments",
                             "exp_results", "Sept_23",
                             paste0("exp_hawkes_param_",
                                    sim_id, ".RDS")))

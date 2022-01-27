#### Exp 12, Jan 26th 2022
#### Investigate whether the initialization function 
#### leads to an improved clustering performance

library(here)

source(here("Experiments/", "utils.R"))
source(here("functions/init_fcn.R"))


### simulate some data, see how the performance changes with the
### init function


Time <- 200
no_sims <- 20
dT_vec <- 1
inter_T <- 1
K_vec <- 2:6
m_vec <- rep(c(100, 200, 400), 2)
sparsity <- 0.8 # prop of edges which can have events

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

model <- "Hawkes"
if(sim_id <= 3){
  model = "Poisson"
}

results <- list()
m <- m_vec[sim_id]

for(exp_num in seq_along(K_vec)) {
  K <- K_vec[exp_num]
  dT <- 1
  curr_dt_sims <- tibble()
  cat("Current K:", K, "\n")
  cat("Current m:", m, "\n")
  cat(model, "\n")
  for(sim in 1:no_sims){
    cat("Sim:", sim, "======\n")
    ## baseline rate of the process
    true_Mu <- matrix(0.05, 
                      nrow = K, ncol = K, byrow = T)
    diag(true_Mu) <- c(rep(0.5, K-1), 1)
    ## excitation, if used (for Hawkes)
    true_B <- matrix(0, nrow = K, ncol = K, byrow = TRUE)
    diag(true_B) <- 0.5
    if(model == "Poisson") {
      true_B <- matrix(0, K, K)
    }
    Pi <- matrix(1/K, 1, K)
    
    Z <- sample(0:(K-1), size = m, prob = Pi, replace = TRUE)
    # then generate A
    A <- list()
    for(i in 1:m){
      # could sample these here with SBM structure...
      num_edge = m * sparsity
      edge <- sample(m, num_edge) - 1
      edge <- sort(edge)
      A[[i]] <- edge
    }
    alltimes <- sampleBlockHak(Time, A, Z, Mu = true_Mu, B = true_B, lam = 1)
    
    ### then estimate the fits here in each case
    if(model == "Poisson") {
      
      ### run init algorithm
      result <- dense_poisson(alltimes, K)
      Mu_est <- result$est_B
      ## need to pass the estimated clustering also
      init_tau <- matrix(0, nrow = m, ncol = K)
      for(i in seq_along(result$est_clust)){
        init_tau[i, result$est_clust[i]] <- 1
      }
      
      ### will need to modify to account for the decreased number
      ### of events also...
      results_online <- estimate_Poisson_init(full_data = result$rest_events,
                                         A,
                                         m,
                                         K,
                                         Time,
                                         dT = dT,
                                         B = Mu_est,
                                         inter_T,
                                         init_tau,
                                         start = result$cut_off,
                                         is_elbo = FALSE)
    }
    if(model == "Hawkes"){
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
                                                     B, Mu, tau,
                                                     inter_T, 
                                                     is_elbo = FALSE)
    }
    
    z_est <- apply(results_online$tau, 1, which.max)
    clust_est <- aricode::ARI(Z, z_est)
    ### then save dT, clust_est, m, model
    curr_sim <- tibble(ARI = clust_est,
                       K = K, 
                       nodes = m,
                       model = model)
    curr_dt_sims <- curr_dt_sims %>% 
      bind_rows(curr_sim)
  }
  results[[exp_num]] <- curr_dt_sims 
}


### then save these somewhere
saveRDS(results, file = here("Experiments",
                             "exp_results",
                             paste0("exp_12_", sim_id, ".RDS")))

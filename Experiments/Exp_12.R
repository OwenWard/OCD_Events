#### Exp 12, March 26th 2022
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
inter_T <- 1
K <- 2
m_vec <- c(rep(100, 5), rep(200, 5), rep(400, 5))
sparsity <- 0.5 # prop of edges which can have events

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

model <- "Hawkes"
if(sim_id <= 3){
  model = "Poisson"
}

results <- list()
m <- m_vec[sim_id]
m0_vec <- c( 100*c(1/10, 1/5, 1/4, 1/2, 1),
             200*c(1/10, 1/5, 1/4, 1/2, 1),
             400*c(1/10, 1/5, 1/4, 1/2, 1))
m0_curr <- m0_vec[sim_id]

n0_vals <- seq(from = 5, to = 50, by = 5)

# for(exp_num in seq_along(m_vec)) {
#   dT <- 1
curr_dt_sims <- tibble()
cat("Current K:", K, "\n")
cat("Current m:", m, "\n")
cat(model, "\n")
for(sim in 1:no_sims){
  cat("Sim:", sim, "======\n")
  ## baseline rate of the process
  true_Mu <- matrix(0.05, 
                    nrow = K, ncol = K, byrow = T)
  diag(true_Mu) <- 0.5:K
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
    edge <- sort(edge[edge!= i-1]) ## remove self edges
    A[[i]] <- edge
  }
  alltimes <- sampleBlockHak(Time, A, Z, Mu = true_Mu, B = true_B, lam = 1)
  ### then estimate the fits here in each case
  if(model == "Poisson") {
    
    for(curr_n0 in n0_vals){
      ### run init algorithm
      # result <- dense_poisson(alltimes, K, n0 = curr_n0, m)
      # m0_curr <- m
      result <- sparse_poisson(alltimes, K, n0 = curr_n0, m, m0 = m0_curr)
      # while(sum(is.nan(result$est_B)) > 0) {
      #   result <- dense_poisson(alltimes, K, n0 = curr_n0)
      #   ## just run again to avoid this issue
      # }
      Mu_est <- result$est_B
      ## need to pass the estimated clustering also
      init_tau <- matrix(0, nrow = m, ncol = K)
      for(i in seq_along(result$est_clust)){
        init_tau[i, result$est_clust[i]] <- 1
      }
      ### check the initial ARI
      aricode::ARI(result$est_clust, Z)
      
      ### will need to modify to account for the decreased number
      ### of events also...
      results_online_init <- estimate_Poisson_init(full_data = 
                                                     result$rest_events,
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
      ### compare to not using init function
      B <- matrix(runif(K * K), K, K)
      norm_online <- estimate_Poisson(full_data = alltimes,
                                      A = A,
                                      m,
                                      K,
                                      Time,
                                      dT = 1,
                                      B,
                                      inter_T = 1,
                                      is_elbo = FALSE)
      stan_est <- apply(norm_online$tau, 1, which.max)
      
      z_est <- apply(results_online_init$tau, 1, which.max)
      clust_est_init <- aricode::ARI(Z, z_est)
      
      clust_est_norm <- aricode::ARI(stan_est, Z)
      
      ### then save dT, clust_est, m, model
      curr_sim <- tibble(ARI = c(clust_est_init, clust_est_norm),
                         K = K, 
                         nodes = m,
                         model = model,
                         init = c("Init", "No Init"),
                         n0 = curr_n0,
                         m0 = m0_curr,
                         sparsity = sparsity)
      curr_dt_sims <- curr_dt_sims %>% 
        bind_rows(curr_sim)
    }
  }
}

results <- curr_dt_sims
  
### then save these somewhere
saveRDS(results, file = here("Experiments",
                             "exp_results",
                             paste0("exp_12_m0_sim", sim_id, ".RDS")))

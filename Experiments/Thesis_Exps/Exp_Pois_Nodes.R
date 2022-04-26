#### Exp varying the number of Nodes, April 25th 2022
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
m_vec <- c(rep(100, 3), rep(200, 3), rep(500, 3), rep(1000, 3))

sparsity <- 0.1 # prop of edges which can have events

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

model <- "Poisson"
# if(sim_id <= 3){
#   model = "Poisson"
# }

results <- list()
m <- m_vec[sim_id]
m0_vec <- c( 100*c(1/10, 1/4, 1/2),
             200*c(1/10, 1/4, 1/2),
             500*c(1/10, 1/4, 1/2),
             1000*c(1/10, 1/4, 1/2))
m0_curr <- m0_vec[sim_id]

n0_vals <- c(10, 25, 50)

# for(exp_num in seq_along(m_vec)) {
#   dT <- 1
curr_dt_sims <- tibble()
cat("Current K:", K, "\n")
cat("Current m:", m, "\n")
cat(model, "\n")
for(sim in 1:no_sims){
  cat("Sim:", sim, "======\n")
  ## baseline rate of the process
  true_Mu <- matrix(c(2, 0.05, 0.15, 1.5), 
                    nrow = K, ncol = K, byrow = T)
  # diag(true_Mu) <- c(2, 4)
  # true_Mu[1, 2] <- 0.25
  # true_Mu[2, 2] <- 0.75
  ## excitation, if used (for Hawkes)
  true_B <- matrix(0, nrow = K, ncol = K, byrow = TRUE)
  diag(true_B) <- 0.5
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
      # aricode::ARI(result$est_clust, Z)
      
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
      z_est <- apply(results_online_init$tau, 1, which.max)
      clust_est_init <- aricode::ARI(Z, z_est)
      # print("Init Worked")
      ### then save dT, clust_est, m, model
      curr_sim <- tibble(ARI = clust_est_init,
                         K = K, 
                         nodes = m,
                         model = model,
                         init = "Init",
                         n0 = curr_n0,
                         m0 = m0_curr,
                         sparsity = sparsity)
      curr_dt_sims <- curr_dt_sims %>% 
        bind_rows(curr_sim)
    }
  }
  ### then do random init down here, bind it to curr_dt_sims
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
  (clust_est_norm <- aricode::ARI(stan_est, Z))
  # print("Random Worked")
  curr_sim_rand <- tibble(ARI = clust_est_norm,
                          K = K, 
                          nodes = m,
                          model = model,
                          init = "No Init",
                          n0 = NA,
                          m0 = NA,
                          sparsity = sparsity)
  curr_dt_sims <- curr_dt_sims %>% 
    bind_rows(curr_sim_rand)
  
}

results <- curr_dt_sims

### then save these somewhere
saveRDS(results, file = here("Experiments",
                             "thesis_output",
                             paste0("exp_pois_nodes_april_25_rho_",
                                    100*sparsity, sim_id, ".RDS")))

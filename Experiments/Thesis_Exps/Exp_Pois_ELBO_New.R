#### Exp varying the number of Nodes, April 25th 2022
#### Investigate whether the initialization function 
#### leads to an improved clustering performance
.libPaths("/moto/stats/users/ogw2103/rpackages")
library(here)

source(here("Experiments/", "utils.R"))
source(here("functions/init_fcn.R"))


### simulate some data, see how the performance changes with the
### init function


no_sims <- 50
dT <- 1
inter_T <- 1
K <- 2
# m_vec <- c(rep(100, 3), rep(200, 3), rep(500, 3), rep(1000, 3),
#            rep(5000, 3))
Time <- 200
m_vec <- c(50, 100, 200, 500, 1000)

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


m0_curr <- m/4
n0_vals <- 5

# for(exp_num in seq_along(m_vec)) {
#   dT <- 1
cat("Current K:", K, "\n")
cat("Current m:", m, "\n")
cat(model, "\n")
for(sim in 1:no_sims){
  cat("Sim:", sim, "======\n")
  ## baseline rate of the process
  true_Mu <- matrix(c(2, 0.25, 0.05, 1.5), 
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
  Pi <- c(0.5, 0.5)
  
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
                                                   is_elbo = TRUE)
      z_est <- apply(results_online_init$tau, 1, which.max)
      clust_est_init <- aricode::ARI(Z, z_est)
      # print("Init Worked")
      ### then save dT, clust_est, m, model
      curr_sim <- list(ARI = clust_est_init,
                       K = K, 
                       nodes = m,
                       model = model,
                       init = "Init",
                       n0 = curr_n0,
                       m0 = m0_curr,
                       sparsity = sparsity,
                       Time = Time,
                       elbo = results_online_init$full_ELBO)
      results[[sim]] <- curr_sim
    }
  }
}


### then save these somewhere
saveRDS(results, file = here("Experiments",
                             "thesis_output",
                             paste0("exp_pois_elbo_nodes_may_19_rho_",
                                    100*sparsity, sim_id, ".RDS")))

## fixed here meaning fixed n0, m0
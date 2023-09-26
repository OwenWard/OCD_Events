##### Examine compute time as we vary dT



.libPaths("/moto/stats/users/ogw2103/rpackages")
library(here)

source(here("functions/", "utils.R"))
source(here("functions/init_fcn.R"))
library(rbenchmark)

### simulate some data, see how the performance changes with the
### init function


Time <- 200
no_sims <- 50
dT <- 1
inter_T <- 1
K <- 2

m <- 100

sparsity <- 0.15 # prop of edges which can have events

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

dT_vec <- seq(from = 0.25, to = 5, by = 0.25)
## 50 scenarios here now

sparsity <- 0.15
no_sims <- 50
curr_dT <- dT_vec[sim_id]
results <- list()
n0 <- 40
m0 <- m/10


# for(exp_num in seq_along(m_vec)) {
#   dT <- 1
window <- 10
curr_dt_sims <- tibble()
cat("Current K:", K, "\n")
cat("Current m:", m, "\n")
cat("Current window:", window, "\n")

for(sim in 1:no_sims){
  cat("Sim:", sim, "======\n")
  ## baseline rate of the process
  
  ####
  ## simulate from inhomogeneous Poisson instead here
  H <- 2
  MuA <- array(0, c(K, K, H))
  MuA[, , 1] <- matrix(c(0.8, 0.2, 0.6, 0.4), 2, 2)
  MuA[, , 2] <- matrix(c(0.4, 0.7, 0.2, 0.7), 2, 2)
  # 
  # MuA[, , 1] <- matrix(c(0.08, 0.02, 0.06, 0.04)/2, 2, 2)
  # MuA[, , 2] <- matrix(c(0.04, 0.07, 0.02, 0.07)/2, 2, 2)
  ####
  
  ## excitation, if used (for Hawkes)
  true_B <- matrix(0, nrow = K, ncol = K, byrow = TRUE)
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
  alltimes <- sampleBlockHak_nonhomo(Time,
                                     A,
                                     Z,
                                     MuA = MuA,
                                     B = true_B, 
                                     window = window,
                                     lam = 1)
  print("Simulated Data")
  ### then estimate the fits here in each case
  
  
  ### then do random init down here, bind it to curr_dt_sims

  # result <- dense_inhom_Poisson(alltimes, K,
  #                               H = H,
  #                               window = curr_wind,
  #                               t_start = 0,
  #                               n0 = n0, m)
  ## to check with sparse also
  sim_time <- benchmark(sim = {
  result <- sparse_inhom_Poisson(alltimes,
                                 K,
                                 H,
                                 window = window,
                                 t_start = 0,
                                 n0 = n0,
                                 m, m0)
  Mu_est <- result$est_Mu
  ## need to pass the estimated clustering also
  init_tau <- matrix(0, nrow = m, ncol = K)
  for(i in seq_along(result$est_clust)){
    init_tau[i, result$est_clust[i]] <- 1
  }
  ### check the initial ARI
  cat("Init Scheme \n")
  print(aricode::ARI(result$est_clust, Z))
  cat("-------\n")
  #     
  #     ### will need to modify to account for the decreased number
  #     ### of events also...
  results_online_init <- nonhomo_Pois_est_init(
    alltimes = result$rest_events,
    A,
    m,
    K,
    H,
    window,
    Time,
    dT = curr_dT,
    gravity = 0.01,
    MuA_start = Mu_est,
    init_tau,
    start = result$cut_off,
    full_data = result$rest_events, ## won't be used, doesn't matter
    is_elbo = FALSE)}, replications = 1)
  z_est <- apply(results_online_init$tau, 1, which.max)
  clust_est_init <- aricode::ARI(Z, z_est)
  cat("Post Init \n")
  print(clust_est_init)
  print(table(z_est, Z))
  cat("------\n")
  
  # print("Random Worked")
  curr_sim <- tibble(ARI = clust_est_init,
                          K = K, 
                          nodes = m,
                          time = sim_time$elapsed,
                          dT = curr_dT,
                          window = window,
                          sparsity = sparsity,
                          winds = (Time - n0)/curr_dT,
                          sim = sim)
  curr_dt_sims <- curr_dt_sims %>% 
    bind_rows(curr_sim)

}

results <- curr_dt_sims


saveRDS(results, file = here("Experiments",
                             "exp_results", "Sept_23",
                             paste0("exp_inhom_pois_time_",
                                    curr_dT, "_", sim_id, ".RDS")))

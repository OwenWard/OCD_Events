#### Experiment InPois Window, Sept 22nd
#### to investigate the relationship between 
#### the window size and the clustering result,
#### for a fixed data set and initial B matrix

.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)

source(here("functions", "utils.R"))
source(here("functions/init_fcn.R"))

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

dT_vec <- seq(from = 0.25, to = 5, by = 0.25)
## 50 scenarios here now

sparsity <- 0.15
no_sims <- 50
curr_dT <- dT_vec[sim_id]


results <- tibble()

for(sim in 1:no_sims) {
  cat("Sim:", sim, "\n")
  m <- 200
  Time <- 200
  n0 <- 40
  m0 <- m/10
  inter_T <- 1
  K <- 2
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
  ###
  result <- sparse_inhom_Poisson(alltimes,
                                 K,
                                 H,
                                 window = curr_dT,
                                 t_start = 0,
                                 n0 = n0,
                                 m, m0)
  Mu_est <- result$est_Mu
  ## need to pass the estimated clustering also
  init_tau <- matrix(0, nrow = m, ncol = K)
  for(i in seq_along(result$est_clust)){
    init_tau[i, result$est_clust[i]] <- 1
  }
  a <- capture.output(results_online_init <- nonhomo_Pois_est_init(
    alltimes = result$rest_events,
    A,
    m,
    K,
    H,
    curr_dT,
    Time,
    dT = dT,
    gravity = 0.01,
    MuA_start = Mu_est,
    init_tau,
    start = result$cut_off,
    full_data = result$rest_events, ## won't be used, doesn't matter
    is_elbo = FALSE))
  z_est <- apply(results_online_init$tau, 1, which.max)
  clust_est_init <- aricode::ARI(Z, z_est)
  cat("Post Init \n")
  print(clust_est_init)
  print(table(z_est, Z))
  cat("------\n")
  # print("Init Worked")
  ### then save dT, clust_est, m, model
  curr_sim <- tibble(ARI = clust_est_init,
                     K = K,
                     nodes = m,
                     # model = model,
                     init = "Init",
                     n0 = n0,
                     window = curr_dT,
                     # m0 = m0_curr,
                     sparsity = sparsity,
                     sim = sim)
  ###
  results <- results %>%  bind_rows(curr_sim)
}

# all_results <- list(results, init)

saveRDS(results, file = here("Experiments",
                             "exp_results", "Sept_23",
                             paste0("exp_pois_wind_", curr_dT, ".RDS")))
# to save it after each run also


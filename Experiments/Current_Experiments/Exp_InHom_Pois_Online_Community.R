#### Exp verifying online community recovery, Apr 26th
#### Investigate whether the initialization function 
#### leads to an improved clustering performance
.libPaths("/moto/stats/users/ogw2103/rpackages")
library(here)

source(here("functions", "utils.R"))
source(here("functions/init_fcn.R"))


### simulate some data, see how the performance changes with the
### init function


Time <- 200
no_sims <- 50
dT <- 1
inter_T <- 20
K <- 2
m_vec <- c(100, 200, 500, 1000)

sparsity <- 0.15 # prop of edges which can have events

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

results <- list()
m <- m_vec[sim_id]


n0 <- 40
m0 <- m/2

window <- 10
curr_dt_sims <- tibble()
cat("Current K:", K, "\n")
cat("Current m:", m, "\n")
for(sim in 1:no_sims){
  cat("Sim:", sim, "======\n")
  ## baseline rate of the process
  H <- 2
  MuA <- array(0, c(K, K, H))
  MuA[, , 1] <- matrix(c(0.8, 0.2, 0.6, 0.4), 2, 2)
  MuA[, , 2] <- matrix(c(0.4, 0.7, 0.2, 0.7), 2, 2)
  
  true_B <- matrix(0, nrow = K, ncol = K)
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
  for(curr_wind in 1:10) {
    # result <- dense_inhom_Poisson(alltimes, K,
    #                               H = H,
    #                               window = curr_wind,
    #                               t_start = 0,
    #                               n0 = n0, m)
    ## to check with sparse also
    result <- sparse_inhom_Poisson(alltimes,
                                   K,
                                   H,
                                   window = curr_wind,
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
    init_ari <- aricode::ARI(result$est_clust, Z)
    print(init_ari)
    cat("-------\n")
    a <- capture.output(results_online_init <- nonhomoPois_est_init(
      alltimes = result$rest_events,
      A,
      m,
      K,
      H,
      curr_wind,
      Time,
      dT = dT,
      gravity = 0.01,
      MuA_start = Mu_est,
      init_tau,
      start = result$cut_off,
      is_elbo = FALSE))
    ### need to capture all the output along the way here
    keep_inter <- ((result$cut_off+1):Time)%%10 == 0
    keep_inter <- keep_inter[-length(keep_inter)]
    ## drop the last which we already compute below
    inter_times <- ((result$cut_off+1):Time)[keep_inter]
    ##
    inter_tau <- results_online_init$inter_tau[, ,keep_inter]
    inter_z <- apply(inter_tau, c(1, 3), which.max)
    inter_ari <- apply(inter_z, 2, function(x) aricode::ARI(x, Z))
    
    ###
    z_est <- apply(results_online_init$tau, 1, which.max)
    clust_est_final <- aricode::ARI(Z, z_est)
    cat("Post Init \n")
    print(clust_est_final)
    print(table(z_est, Z))
    cat("------\n")
    # print("Init Worked")
    ### then save dT, clust_est, m, model
    curr_sim <- tibble(ARI = c(init_ari, inter_ari, clust_est_final),
                       K = K,
                       time = c(n0, inter_times, Time),
                       init = "Init",
                       n0 = n0,
                       window = curr_wind,
                       m0 = m0,
                       sparsity = sparsity,
                       sim = sim)
    curr_dt_sims <- curr_dt_sims %>%
      bind_rows(curr_sim)
    
  }
}

results <- curr_dt_sims

### then save these somewhere
saveRDS(results, file = here("Experiments",
                             "thesis_output",
                             paste0("exp_in_pois_online_nov_9_rho_",
                                    100*sparsity, sim_id, ".RDS")))

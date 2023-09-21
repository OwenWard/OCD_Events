#### Exp varying the number of Nodes, Oct 10th 2022
#### Investigate performance without initialization
#### for inhomogeneous Poisson first
.libPaths("/moto/stats/users/ogw2103/rpackages")
library(here)

source(here("functions/", "utils.R"))
source(here("functions/init_fcn.R"))


### simulate some data, see how the performance changes with the
### init function


Time <- 200
no_sims <- 50
dT <- 1
inter_T <- 1
K <- 2

m_vec <- c(1000, 5000, 10000, 50000)

sparsity <- 0.01 # prop of edges which can have events

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid


results <- list()
m <- m_vec[sim_id]
# m0_vec <- c( 100*c(1/10, 1/4, 1/2),
#              200*c(1/10, 1/4, 1/2),
#              500*c(1/10, 1/4, 1/2),
#              1000*c(1/10, 1/4, 1/2))
# 
# m0_curr <- m0_vec[sim_id]
# 
# 
# m0_curr <- m/4
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
  for(curr_wind in 1) {
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
    print(aricode::ARI(result$est_clust, Z))
    cat("-------\n")
    #     
    #     ### will need to modify to account for the decreased number
    #     ### of events also...
    a <- capture.output(results_online_init <- nonhomo_Pois_est_init(
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
                       window = curr_wind,
                       # m0 = m0_curr,
                       sparsity = sparsity,
                       sim = sim)
    curr_dt_sims <- curr_dt_sims %>%
      bind_rows(curr_sim)
    MuA_init <- array(runif(K * K * H), dim = c(K, K, H) )
    tau_init <- matrix(1/K, nrow = m, ncol = K)
    
    a <- capture.output(norm_online <- nonhomoPois_estimator(alltimes, 
                                         A,
                                         m,
                                         K,
                                         H, 
                                         window = curr_wind,
                                         T = Time,
                                         dT,
                                         gravity = 0.01,
                                         MuA_start = MuA_init,
                                         tau_start = tau_init))
    
    
    stan_est <- apply(norm_online$tau, 1, which.max)
    (clust_est_norm <- aricode::ARI(stan_est, Z))
    # print("Random Worked")
    curr_sim_rand <- tibble(ARI = clust_est_norm,
                            K = K, 
                            nodes = m,
                            init = "No Init",
                            n0 = NA,
                            window = curr_wind,
                            # m0 = NA,
                            # init_num = sim_id,
                            sparsity = sparsity,
                            sim = sim)
    curr_dt_sims <- curr_dt_sims %>% 
      bind_rows(curr_sim_rand)
  }
}

results <- curr_dt_sims

### then save these somewhere
saveRDS(results, file = here("Experiments",
                             "exp_results", "June_23",
                             paste0("exp_inpois_nodes_fixed_june_22_rho_",
                                    100 * sparsity, "_", sim_id, ".RDS")))

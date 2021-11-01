#### Experiment 2, Oct 26th
#### to investigate the relationship between 
#### the window size and the clustering result,
#### for a fixed data set and initial B matrix

.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)

source(here("Experiments/", "utils.R"))

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

dT_vec <- seq(from = 0.1, to = 5, by = 0.1)
## 50 scenarios here now


no_sims <- 20
curr_dT <- dT_vec[sim_id]


### then load in a dataset and an initial B? one which works 
### well for a known dT
# exp <- readRDS(here("Experiments/", "sim_pars.RDS"))
# 
# B <- exp$init_B
# sim1 <- exp$sim
### now random B and random dataset each time also

results <- list()
init <- list()


for(sim in 1:no_sims) {
  cat("Sim:", rep, "\n")
  n <- 100
  Time <- 100
  intens1 <- c(2)
  intens2 <- c(1)
  true_B <- matrix(c(intens1, 0.05, 0.05, intens2), 
                   nrow = 2, ncol = 2, byrow = T)
  
  intens <- matrix(c(intens1, 0.05, 0.05, intens2), 4, 1)
  sim1 <- generateDynppsbmConst(intens = intens,
                                            Time = Time,
                                            n = n, 
                                            prop.groups = c(0.5, 0.5))
  proc_sim <- format_sims(sim_data = sim1, n = n)
  
  K <- 2
  m <- n
  B <- matrix(runif(K * K), K, K)
  dT <- dT_vec[sim]
  inter_T <- 1
  
  results_online <- estimate_Poisson(full_data = proc_sim$events,
                                     A = proc_sim$edge,
                                     m,
                                     K,
                                     Time,
                                     dT,
                                     B,
                                     inter_T,
                                     is_elbo = TRUE)
  
  # compute rand index
  z_true <- apply(sim1$z, 2, function(x) which(x == 1))
  z_est <- apply(results_online$tau, 1, which.max)
  clust_est <- aricode::ARI(z_true, z_est)
  
  sim_pars <- list(
    B = B, 
    est_elbo = results_online$AveELBO,
    clust = clust_est,
    tau = results_online$tau,
    dT = dT
  )
  results[[sim]] <- sim_pars
}

# all_results <- list(results, init)

saveRDS(results, file = here("Experiments",
                                 "exp_results",
                                 paste0("exp_2_", curr_dT, ".RDS")))
# to save it after each run also


# all_results <- list(results, init)
# saveRDS(all_results, file = here("Experiments",
#                              "exp_results",
#                              paste0(Sys.date,"exp_dt_rep.RDS")))

# 
# #### old code for first run ####
# inter_T <- 1
# results <- list()
# 
# for(sim in 1:no_sims) {
#   dT <- dT_vec[sim]
#   
#   results_online <- estimate_Poisson(full_data = proc_sim$events,
#                                      A = proc_sim$edge,
#                                      m,
#                                      K,
#                                      Time,
#                                      dT,
#                                      B,
#                                      tau,
#                                      Pi,
#                                      S,
#                                      inter_T,
#                                      is_elbo = TRUE)
#   
#   # compute rand index
#   z_true <- apply(sim1$z, 2, function(x) which(x == 1))
#   z_est <- apply(results_online$tau, 1, which.max)
#   clust_est <- aricode::ARI(z_true, z_est)
#   
#   sim_pars <- list(
#     dT = dT, 
#     est_eblo = results_online$AveELBO,
#     clust = clust_est
#   )
#   results[[sim]] <- sim_pars
# }
# 
# 
# saveRDS(results, file = here("Experiments", 
#                              "exp_results",
#                              "exp_dt.RDS"))

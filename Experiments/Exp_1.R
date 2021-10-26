#### Experiments, Oct 26th
#### To investigate the relationship between the 
#### the initial B matrix and the clustering result


.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)

source(here("Experiments/", "utils.R"))

nsims <- 100

results <- list()

for(sim in 1:nsims) {
  Time <- 500
  n <- 100
  intens1 <- c(1.2)
  intens2 <- c(0.75)
  true_B <- matrix(c(intens1, 0.15, 0.3, intens2), 
                   nrow = 2, ncol = 2, byrow = T)
  
  intens <- matrix(c(intens1, 0, 0.1, intens2), 4, 1)
  # this is essentially the K*K matrix stretched out as a single col
  
  system.time(sim1 <- generateDynppsbmConst(intens = intens,
                                            Time = Time,
                                            n = n, 
                                            prop.groups = c(0.5, 0.5)))
  
  proc_sim <- format_sims(sim_data = sim1, n = n)
  
  K <- 2
  m <- n
  Pi <- rep(1/K, K)
  B <- matrix(runif(K * K), K, K)
  # diag(B) <- rnorm(K, mean = 1, sd = 0.1)
  tau <- matrix(1, nrow = m, ncol = K)
  tau <- tau/rowSums(tau)
  S <- matrix(1/K, nrow = m, ncol = K)
  
  dT <- 1
  inter_T <- 1
  
  # capture output to not print out
  out <- capture.output(results_online <- estimate_Poisson(full_data = proc_sim$events,
                                     A = proc_sim$edge,
                                     m,
                                     K,
                                     Time,
                                     dT,
                                     B,
                                     tau,
                                     Pi,
                                     S,
                                     inter_T,
                                     is_elbo = TRUE))
  
  # compute rand index
  z_true <- apply(sim1$z, 2, which.max)
  z_est <- apply(results_online$tau, 1, which.max)
  clust_est <- aricode::ARI(z_true, z_est)
  
  sim_pars <- list(
    B = B, 
    est_eblo = results_online$AveELBO,
    clust = clust_est
  )
  results[[sim]] <- sim_pars
}

saveRDS(results, file = here("Experiments",
                             "exp_results",
                             "exp1_long.RDS"))

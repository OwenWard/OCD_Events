#### Poisson Simulation
.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)
source(here("Experiments/", "utils.R"))
# source(here("functions/init_fcn.R"))
## Compare log likelihood, and time

no_sims <- 50

sparsity <- 0.1
Time <- 100
dT <- 1
K <- 3
Mu <- matrix(c(0.6, 0.2, 0.3, 0.1, 1.0, 0.4, 0.5, 0.4, 0.8), K, K, byrow = TRUE)
B <- matrix(0, K, K, byrow = TRUE)
diag(B) <- 0.5
m <- 100
Pi <- matrix(c(0.4, 0.3, 0.3), 1, 3)
Z <- c(rep(0, m * Pi[1]), rep(1, m * Pi[2]), rep(2, m * Pi[3]))

A <- list()
for(i in 1:m){
  # could sample these here with SBM structure...
  num_edge = m * sparsity
  edge <- sample(m, num_edge) - 1
  edge <- sort(edge[edge!= i-1]) ## remove self edges
  A[[i]] <- edge
}

system.time(alltimes <- sampleBlockHak(T, A, Z, Mu, B, lam = 1))

results <- tibble()

for(sim in 1:no_sims){
  cat("Sim:", sim, "-----\n")
  init_tau <- matrix(1/K, nrow = m, ncol = K)
  init_B <- matrix(runif(K * K), nrow = K, ncol = K)
  init_Mu <- matrix(runif(K * K), nrow = K, ncol = K)
  
  
  t1 <- bench::mark(
    results_hawkes <- online_estimator_eff_revised(alltimes, 
                                                       A,
                                                       m,
                                                       K,
                                                       T = Time,
                                                       dT, 
                                                       lam = 0.5,
                                                       init_B,
                                                       init_Mu,
                                                       init_tau,
                                                       inter_T = 1),
    iterations = 1)
  
  ### then compute the batch estimator
  t2 <- bench::mark(
    results_hawkes_batch <- batch_estimator(alltimes, 
                                            A,
                                            m,
                                            K,
                                            T = Time,
                                            dT,
                                            lam = 0.5,
                                            init_B,
                                            init_Mu,
                                            init_tau,
                                            itermax = 100,
                                            stop_eps = 0.01),
    iterations = 1)
  ### compare the log likelihood of each
  Z_online <- apply(results_hawkes$tau, 1, which.max) - 1
  ll_online <- get_loglik_Hak(alltimes,
                              t_start = 0, t_end = Time,
                              Z = Z_online, Mu = results_hawkes$Mu,
                              B = results_hawkes$B, Pi = results_hawkes$Pi,
                              A, results_hawkes$lam, m, K)
  Z_batch <- apply(results_hawkes_batch$tau, 1, which.max) - 1
  ll_batch <- get_loglik_Hak(alltimes, 
                             t_start = 0, t_end = Time,
                             Z = Z_batch, Mu = results_hawkes_batch$Mu,
                             B = results_hawkes_batch$B, 
                             Pi = results_hawkes_batch$Pi,
                             A, results_hawkes_batch$lam, m, K)
  
  ### store everything here
  diff = 1 - abs(ll_batch - ll_online)/abs(ll_batch)
  
  curr_sim <- tibble(
    sim = sim,
    t_on = t1$total_time,
    t_batch = t2$total_time,
    ll_online = ll_online,
    ll_batch = ll_batch,
    perc_diff = diff
  )
  results <- results %>% 
    bind_rows(curr_sim)
}



## then save results to something down here
saveRDS(results, file = here("Experiments",
                             "thesis_output",
                             paste0("exp_time_hawkes", ".RDS")))

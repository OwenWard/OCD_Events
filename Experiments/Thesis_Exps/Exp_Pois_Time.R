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

system.time(alltimes <- sampleBlockHak(Time, A, Z, Mu, B, lam = 1))

result_sims <- tibble()

for(sim in 1:no_sims){
  cat("Sim:", sim, "-----\n")
  init_B <- matrix(runif(K * K), nrow = K, ncol = K)
  
  
  t1 <- bench::mark(results <- estimate_Poisson(full_data = alltimes,
                                                A, m, K, Time, dT,
                                                init_B, inter_T = 1),
                    iterations = 1)
  
  ### then compute the batch estimator
  t2 <- bench::mark(results_pois_batch <- batch_estimator_hom_Poisson(
    alltimes,
    A,
    m,
    K,
    Time,
    itermax = 100,
    stop_eps = 0.01),
  iterations = 1)
  ### compare the log likelihood of each
  ll_online <- computeLL(alltimes, results$tau, results$B,
                         results$Pi, A, m, K, currT = Time)
  ll_batch <- computeLL(alltimes, results_pois_batch$tau,
                        results_pois_batch$B, results_pois_batch$Pi,
                        A, m , K, currT = Time)
  
  ### store everything here
  diff = 1 - abs(ll_batch - ll_online)/abs(ll_batch)
  
  curr_sim <- tibble(
    sim = sim,
    t_on = as.numeric(t1$total_time),
    t_batch = as.numeric(t2$total_time),
    ll_online = ll_online,
    ll_batch = ll_batch,
    perc_diff = diff
  )
  result_sims <- result_sims %>% 
    bind_rows(curr_sim)
}



## then save results to something down here
saveRDS(results, file = here("Experiments",
                             "thesis_output",
                             paste0("exp_time_pois", ".RDS")))

#### Experiments, Nov 1st
#### To explore the relationship between result and the predicted loss
#### in terms of expected events in the next time window
.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)

source(here("Experiments/", "utils.R"))

nsims <- 10

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid
Times <- c(50, 100, 200, 500)
Time <- Times[sim_id]

results <- list()

for(sim in 1:nsims) {
  cat("Sim:", sim, "\n")
  # Time <- 100
  n <- 100
  intens1 <- c(2)
  intens2 <- c(1)
  intens <- matrix(c(intens1, 0.05, 0.05, intens2), 4, 1)
  true_B <- matrix(c(intens1, 0.05, 0.05, intens2), 
                   nrow = 2, ncol = 2, byrow = T)
  # this is essentially the K*K matrix stretched out as a single col
  system.time(sim1 <- generateDynppsbmConst(intens = intens,
                                            Time = Time,
                                            n = n, 
                                            prop.groups = c(0.5, 0.5)))
  
  proc_sim <- format_sims(sim_data = sim1, n = n)
  K <- 2
  m <- n
  B <- matrix(runif(K * K), K, K)
  dT <- 1
  inter_T <- 1
  # capture output to not print out
  results_online <- estimate_Poisson(full_data = proc_sim$events,
                                     A = proc_sim$edge,
                                     m,
                                     K,
                                     Time,
                                     dT,
                                     B,
                                     inter_T,
                                     is_elbo = TRUE)
  
  
  B_ests <- results_online$inter_B
  z_true <- apply(sim1$z, 2, which.max)
  out <- compute_regret(full_data = proc_sim$events,
                        A = proc_sim$edge, 
                        m,
                        K,
                        Time,
                        dT,
                        true_z = z_true,
                        B_ests = B_ests,
                        tau_ests = results_online$early_tau,
                        true_B = true_B)
  
  
  card_A <- as_tibble(proc_sim$events) %>% 
    select(V1, V2) %>% 
    distinct() %>% 
    nrow()
  ## is this regret function correct?
  est_loss <- -out$EstLLH/card_A
  best_loss <- -out$TrueLLH/card_A
  regret <- cumsum(est_loss) - cumsum(best_loss)

  ### compute batch estimates and online loss estimates here
  onl_loss <- tibble(cumsum(out$Online_Loss_True)/(1:Time),
              cumsum(out$Online_Loss_Est)/(1:Time),
              dT = 1:Time) 
  colnames(onl_loss) <- c("True_Z", "Est_Z", "dT")
  
  tidy_loss <- onl_loss %>% 
    pivot_longer(cols = c(True_Z:Est_Z),
                 values_to = "Loss",
                 names_to = "Z") 
  
  
  
  Dmax <- 2^3
  Nijk <- statistics(sim1$data,
                     n,
                     Dmax,
                     directed=TRUE)
  sol.kernel <- mainVEM(list(Nijk = Nijk, Time = Time),
                        n = m,
                        Qmin  = 2,
                        Qmax = 2,
                        directed = TRUE,
                        method = 'hist',
                        d_part = 0,
                        n_perturb = 0)
  
  b <- exp(sol.kernel[[1]]$logintensities.ql[, 1])
  # these largely match the true intensities
  # sol.kernel[[1]]$tau
  z_vem <- apply(sol.kernel[[1]]$tau, 2, which.max)
  
  batch_b <- matrix(b, nrow = 2, ncol = 2, byrow = 2)
  
  vem_loss <- batch_loss(full_data = proc_sim$events,
                         A = proc_sim$edge, m,
                         K,
                         Time,
                         dT, 
                         batch_B = batch_b, true_z = z_vem)
  batch_average <- mean(vem_loss$Batch_loss/1:Time)
  #####
  z_est <- apply(results_online$tau, 1, which.max)
  clust_est <- aricode::ARI(z_true, z_est)
  
  ave_pred_ll <- tibble(online = as.vector(out$Ave_Pred_LL),
                        batch = as.vector(vem_loss$Batch_Ave_Pred_LL)) %>% 
    mutate(dT = row_number()) %>%
    pivot_longer(cols = online:batch,
                 names_to = "method",
                 values_to = "loglik")
  
  sim_pars <- list(
    B = B,
    z_true = z_true,
    z_est = z_est,
    B_ests = B_ests,
    est_elbo = results_online$AveELBO,
    clust = clust_est,
    regret = regret,
    card_A = card_A,
    batch_ave_loss = batch_average,
    online_loss = tidy_loss,
    pred_llh = ave_pred_ll,
    Time = Time
  )
  results[[sim]] <- sim_pars
}

saveRDS(results, file = here("Experiments",
                             "exp_results",
                             paste0("exp4_", Time, ".RDS")))


### create some nicer plots for these regret functions
### for a given dataset, draw new B and plot the loss each time.


tidy_loss %>% 
  ggplot(aes(dT, Loss, colour = Z)) +
  geom_line() +
  geom_hline(aes(yintercept = batch_average, linetype = "Batch Average")) +
  labs(linetype = "", colour = "Labels Used")

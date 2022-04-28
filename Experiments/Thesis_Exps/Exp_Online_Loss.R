#### Experiments, Nov 1st
#### To explore the relationship between result and the predicted loss
#### in terms of expected events in the next time window
.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)
library(ppsbm)
source(here("Experiments/", "utils.R"))
source(here("functions/init_fcn.R"))

nsims <- 25

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
  # B <- matrix(runif(K * K), K, K)
  dT <- 1
  inter_T <- 1
  # capture output to not print out
  ### use initialization scheme here first also
  result <- dense_poisson(alltimes = proc_sim$events,
                          K = K,
                          n0 = 5,
                          m = m)
  Mu_est <- result$est_B
  init_tau <- matrix(0, nrow = m, ncol = K)
  for(i in seq_along(result$est_clust)){
    init_tau[i, result$est_clust[i]] <- 1
  }
  
  ## modify this below then also
  results_online_init <- estimate_Poisson_init(full_data = 
                                                 result$rest_events,
                                               A = proc_sim$edge,
                                               m,
                                               K,
                                               Time,
                                               dT = dT,
                                               B = Mu_est,
                                               inter_T,
                                               init_tau,
                                               start = result$cut_off,
                                               is_elbo = TRUE)
  ###
  # results_online <- estimate_Poisson(full_data = proc_sim$events,
  #                                    A = proc_sim$edge,
  #                                    m,
  #                                    K,
  #                                    Time,
  #                                    dT,
  #                                    B,
  #                                    inter_T,
  #                                    is_elbo = TRUE)
  
  B_ests <- results_online_init$inter_B
  z_true <- apply(sim1$z, 2, which.max)
  ## need to modify the regret function also
  init_time = Time - result$cut_off
  out <- compute_regret(full_data = 
                          result$rest_events,
                        A = proc_sim$edge, 
                        m,
                        K,
                        init_time,
                        dT,
                        true_z = z_true,
                        B_ests = B_ests,
                        tau_ests = results_online_init$early_tau,
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
  onl_loss <- tibble(cumsum(out$Online_Loss_True)/((Time - init_time + 1):Time),
                     cumsum(out$Online_Loss_Est)/((Time - init_time + 1):Time),
                     dT = (Time - init_time + 1):Time) 
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
  
  batch_b <- matrix(b, nrow = 2, ncol = 2, byrow = TRUE)
  
  vem_loss <- batch_loss(full_data = proc_sim$events,
                         A = proc_sim$edge, m,
                         K,
                         Time,
                         dT, 
                         batch_B = batch_b, true_z = z_vem)
  batch_average <- mean(vem_loss$Batch_loss/1:Time)
  #####
  z_est <- apply(results_online_init$tau, 1, which.max)
  clust_est <- aricode::ARI(z_true, z_est)
  
  batch_pred_ll <- as.vector(vem_loss$Batch_Ave_Pred_LL)[(Time - 
                                                            init_time+1):Time]
  
  ave_pred_ll <- tibble(online = as.vector(out$Ave_Pred_LL),
                        ### need to remove some at the start here
                        ### to compare the batch
                        batch = as.vector(na.omit(batch_pred_ll))) %>% 
    mutate(dT = row_number()) %>%
    pivot_longer(cols = online:batch,
                 names_to = "method",
                 values_to = "loglik")
  
  sim_pars <- list(
    B = Mu_est,
    z_true = z_true,
    z_est = z_est,
    B_ests = B_ests,
    est_elbo = results_online_init$AveELBO,
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
                             "thesis_output",
                             paste0("exp_pois_online_loss_april_28_",
                                    Time, ".RDS")))


### create some nicer plots for these regret functions
### for a given dataset, draw new B and plot the loss each time.


# tidy_loss %>% 
#   ggplot(aes(dT, Loss, colour = Z)) +
#   geom_line() +
#   geom_hline(aes(yintercept = batch_average, linetype = "Batch Average")) +
#   labs(linetype = "", colour = "Labels Used")

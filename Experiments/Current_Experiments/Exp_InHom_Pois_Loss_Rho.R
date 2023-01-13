#### Experiments, Nov 1st
#### To explore the relationship between result and the predicted loss
#### in terms of expected events in the next time window
.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)
library(ppsbm)
source(here("functions", "utils.R"))
source(here("functions", "init_fcn.R"))

nsims <- 100

# prob_edge <- 0.5

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid
# Times <- c(50, 100, 200, 500)
Time <- 200
# Time <- Times[sim_id]

sparsity <- c(0.1, 0.25, 0.5, 0.75, 1)

rho <- sparsity[sim_id]

if(Time == 500){
  nsims <- 10
}

results <- list()

### modify to just use a single simulated dataset 
### and so common "batch optimum"
n <- 100

H <- 2

intens1 <- c(2)
intens2 <- c(1)

intens_c1 <- c(intens1, 0.05, 0.05, intens2)
intens_c2 <- c(intens1/2, 0.05, 0.05, intens2/2)

intens <- cbind(intens_c1, intens_c2, intens_c1, intens_c2)
intens <- matrix(cbind(intens_c1, intens_c2), nrow = 4, ncol = Time/10)

window <- Time/ncol(intens)

# intens <- matrix(c(intens1, 0.05, 0.05, intens2), 4, 1)

### make this varying over time, repeat or something

# true_B <- matrix(c(intens1, 0.05, 0.05, intens2), 
#                  nrow = 2, ncol = 2, byrow = T)
# this is essentially the K*K matrix stretched out as a single col
system.time(sim1 <- gen_ppsbm(intens = intens,
                              Time = Time,
                              n = n, 
                              prop.groups = c(0.5, 0.5), 
                              prob_edge = rho))
proc_sim <- format_sims(sim_data = sim1, n = n)

Dmax <- 2^3
Nijk <- statistics(sim1$data,
                   n,
                   Dmax,
                   directed=TRUE)

m <- n
K <- 2
dT <- 1

sol.kernel <- mainVEM(list(Nijk = Nijk, Time = Time),
                      n = m,
                      Qmin  = 2,
                      Qmax = 2,
                      directed = TRUE,
                      method = 'hist',
                      d_part = 0,
                      n_perturb = 0)

b <- exp(sol.kernel[[1]]$logintensities.ql)
## assuming that H = 2, starts and finishes in different bases
b_used <- b[, c(1, ncol(b))]

b_est <- array(as.vector(b_used), dim = c(K, K, H))

# these largely match the true intensities
# sol.kernel[[1]]$tau
z_vem <- apply(sol.kernel[[1]]$tau, 2, which.max)

batch_b <- b_est

vem_loss <- inhom_batch_loss(full_data = proc_sim$events,
                             A = proc_sim$edge, m,
                             K,
                             Time,
                             dT, 
                             batch_B = batch_b,
                             true_z = z_vem, window, H)
batch_average <- mean(vem_loss$Batch_loss/1:Time)

## sum of the losses divided by the number of windows, which 
## is equal to the time here
actual_average <- sum(vem_loss$Batch_loss)/Time

new_batch <- cumsum(vem_loss$Batch_loss)/1:Time

for(sim in 1:nsims) {
  cat("Sim:", sim, "\n")
  # Time <- 100
  
  inter_T <- 1
  
  result <- sparse_inhom_Poisson(alltimes = proc_sim$events,
                                 K = K,
                                 H = H, 
                                 window = window,
                                 t_start = 0, 
                                 n0 = 20,
                                 m,
                                 m0 = m/2)
  Mu_est <- result$est_Mu
  init_tau <- matrix(0, nrow = m, ncol = K)
  for(i in seq_along(result$est_clust)){
    init_tau[i, result$est_clust[i]] <- 1
  }
  
  
  results_online_init <- nonhomo_Pois_est_init(alltimes = result$rest_events,
                                               A = proc_sim$edge, 
                                               m,
                                               K,
                                               H,
                                               window,
                                               Time,
                                               dT = dT,
                                               gravity = 0.01,
                                               MuA_start = Mu_est,
                                               tau_init = init_tau,
                                               start = result$cut_off, 
                                               full_data = proc_sim$events,
                                               is_elbo = TRUE)
  
  
  
  
  
  Mu_ests <- results_online_init$inter_MuA
  true_MuA <- true_MuA <- array(NA, dim = c(K, K, H))
  true_MuA[, , 1] <- matrix(intens[, 1], K, K, byrow = T)
  true_MuA[, , 2] <- matrix(intens[, 2], K, K, byrow = T)
  z_true <- apply(sim1$z, 2, which.max)
  init_time = Time - result$cut_off
  ## need to modify the regret function also
  out <- compute_regret_inhom(full_data = 
                                result$rest_events,
                              A = proc_sim$edge, 
                              m,
                              K,
                              H, 
                              init_time,
                              dT,
                              window,
                              true_z = z_true,
                              MuA_ests = Mu_ests,
                              tau_ests = results_online_init$inter_tau,
                              true_MuA = true_MuA)
  
  colnames(proc_sim$events) <- c("V1", "V2", "V3")
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
    Mu_est = Mu_est,
    z_true = z_true,
    z_est = z_est,
    Mu_ests = Mu_ests,
    est_elbo = results_online_init$elbo,
    clust = clust_est,
    regret = regret,
    card_A = card_A,
    batch_ave_loss = actual_average,
    new_batch_loss = new_batch,
    online_loss = tidy_loss,
    pred_llh = ave_pred_ll,
    Time = Time,
    sparsity = rho
  )
  results[[sim]] <- sim_pars
}



# saveRDS(results, file = here("Experiments",
#                              "exp_results", "November",
#                              paste0("exp_in_pois_online_loss_nov_22_",
#                                     Time, ".RDS")))

if(Time == 500){
  saveRDS(results, file = here("Experiments",
                               "exp_results", "November",
                               paste0("exp_in_pois_online_loss_jan_12_",
                                      Time, "_", sim_id, ".RDS")))
}else{
  saveRDS(results, file = here("Experiments",
                               "exp_results", "November",
                               paste0("exp_in_pois_online_loss_rho_",
                                      sim_id, ".RDS")))
}


### create some nicer plots for these regret functions
### for a given dataset, draw new B and plot the loss each time.


# tidy_loss %>% 
#   ggplot(aes(dT, Loss, colour = Z)) +
#   geom_line() +
#   geom_hline(aes(yintercept = batch_average, linetype = "Batch Average")) +
#   labs(linetype = "", colour = "Labels Used")

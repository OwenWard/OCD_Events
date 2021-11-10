library(here)
source(here("Experiments/", "utils.R"))


### then do a simulation
### want to modify this to look at the impact of dT for a single
### fixed initialization

### save the data, initial values which give good clustering
# exp <- list(sim = sim1, init_B = B)
# saveRDS(exp, file = here("Experiments/", "sim_pars.RDS"))
###


Time <- 100
n <- 100
intens1 <- c(1.2)
intens2 <- c(0.75)
true_B <- matrix(c(intens1, 0, 0.1, intens2), 
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

results_online <- estimate_Poisson(full_data = proc_sim$events,
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
                                   is_elbo = TRUE)

# compute rand index
z_true <- apply(sim1$z, 2, function(x) which(x == 1))
z_est <- apply(results_online$tau, 1, which.max)
aricode::ARI(z_true, z_est)


### repeat for batch model
results_batch <- batch_estimator_hom_Poisson(alltimes = proc_sim$events,
                                             A = proc_sim$edge,
                                             m,
                                             K,
                                             Time,
                                             B,
                                             tau,
                                             itermax = 20,
                                             stop_eps = 0.01)

z_est_batch <- apply(results_batch$tau, 1, which.max)
aricode::ARI(z_true, z_est_batch)


results_batch$ELBO
results_batch$AveELBO

### compare to estimates from fitting ppsbm model,
### compare ELBO etc
### compute regret using this estimate also

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
sol.kernel[[1]]$tau
z_vem <- apply(sol.kernel[[1]]$tau, 2, which.max)

aricode::ARI(z_true, z_vem)


### compute window loss in this setting for batch estimates...
batch_b <- matrix(b, nrow = 2, ncol = 2, byrow = 2)

vem_loss <- batch_loss(full_data = proc_sim$events, A = proc_sim$edge, m,K,Time,dT, 
           batch_B = batch_b, true_z = z_vem)

# is this what I want?
mean(vem_loss$Batch_loss/(1:100))
## this doesn't look correct?
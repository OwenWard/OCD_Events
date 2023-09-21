##### Design Events Used Simulation for Inhom Poisson ####
library(here)
source(here("functions", "utils.R"))
source(here("functions", "init_fcn.R"))

set.seed(100)

m <- 200
Time <- 200
dT <- 1
K <- 2
sparsity <- 0.15

window <- 10

H <- 2
MuA <- array(0, c(K, K, H))
MuA[, , 1] <- matrix(c(0.8, 0.2, 0.6, 0.4), 2, 2)
MuA[, , 2] <- matrix(c(0.4, 0.7, 0.2, 0.7), 2, 2)

true_B <- matrix(0, nrow = K, ncol = K, byrow = TRUE)
Pi <- c(0.5, 0.5)

Z <- sample(0:(K-1), size = m, prob = Pi, replace = TRUE)

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

### then look at this process
n0 <- 20
m0 <- 50

results <- sparse_inhom_Poisson(alltimes, K, H, window, t_start = 0,
                                n0, m, m0)

aricode::ARI(results$est_clust, Z)

Mu_A_est <- results$est_Mu
## need to pass the estimated clustering also
init_tau <- matrix(0, nrow = m, ncol = K)
for(i in seq_along(results$est_clust)){
  init_tau[i, results$est_clust[i]] <- 1
}
### check the initial ARI
aricode::ARI(results$est_clust, Z)


results_online_init <- nonhomo_Pois_est_init(alltimes = results$rest_events,
                                            A, 
                                            m,
                                            K, 
                                            H,
                                            window,
                                            Time,
                                            dT,
                                            gravity = 0.01,
                                            MuA_start = Mu_A_est,
                                            tau_init = init_tau,
                                            start = n0, 
                                            full_data = alltimes,
                                            is_elbo = TRUE)

### what type of elbo is this, what should it look like, compared to the 
## homogeneous Poisson version?
results_online_init$elbo



plot(results_online_init$elbo, type = "l")
table(Z, apply(results_online_init$tau, 1, which.max))

MuA_start <- array(runif(K * K * H), dim = c(K, K, H))
tau_start <- matrix(runif(K * m), nrow = m, ncol = K)

tau_start <- tau_start / rowSums(tau_start)

results_batch <- batch_nonhomoPois_estimator(alltimes,
                                             A,
                                             m,
                                             K,
                                             H,
                                             window,
                                             Time,
                                             dT,
                                             gravity = 0.1,
                                             MuA_start = MuA_start,
                                             tau_start = tau_start,
                                             itermax = 100, stop_eps = 0.01)

aricode::ARI(apply(results_batch$tau, 1, which.max), Z)
## does ok in this example


plot(results_batch$ELBO[results_batch$ELBO!= 0], type = "l")

### this looks fine and useable
init_events <- alltimes[alltimes[,3]< n0, ]

batch_data <- tibble(ELBO = 
                       as.numeric(results_batch$ELBO[1:5]),
                     method = "BATCH") %>% 
  mutate(index = row_number(), events_seen = nrow(alltimes) * index)


online_data <- tibble(ELBO = as.numeric(results_online_init$elbo),
                      method = "ONLINE",
                      events_seen = as.numeric(results_online_init$cum_events) +
                        nrow(init_events)) 



##

axis_title <- 16
axis_text_size <- 14
legend_text <- 14
line_size <- 1

events_plot <- bind_rows(batch_data, online_data) %>% 
  ggplot(aes(events_seen, ELBO, colour = method)) +
  geom_line() +
  scale_x_log10() +
  labs(x = "Events Seen", colour = "Method") +
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = axis_text_size),
        axis.title = element_text(size = axis_title),
        strip.text.x = element_text(size = axis_title),
        legend.text = element_text(size = legend_text)) 

events_plot

ggsave("elbo_events.png",
       events_plot,
       path = here("output_scripts/rev_figure/"),
       dev = "png",
       height = 5,
       width = 7, dpi = 300)


## make this plot in terms of events used instead

new_events_plot <- bind_rows(batch_data, online_data) %>% 
  mutate(perc_events = events_seen/nrow(alltimes)) %>% 
  ggplot(aes(perc_events, ELBO, colour = method)) +
  geom_line() +
  scale_x_continuous(labels = scales::percent, 
                     breaks = scales::pretty_breaks(n = 3)) +
  labs(x = "Events Seen", colour = "Method") +
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = axis_text_size),
        axis.title = element_text(size = axis_title),
        strip.text.x = element_text(size = axis_title),
        legend.text = element_text(size = legend_text)) 

ggsave("elbo_events_2.png",
       new_events_plot,
       path = here("output_scripts/rev_figure/"),
       dev = "png",
       height = 5,
       width = 7, dpi = 300)


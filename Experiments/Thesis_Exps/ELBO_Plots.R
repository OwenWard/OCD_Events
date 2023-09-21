library(here)
source(here("functions", "utils.R"))
source(here("functions","init_fcn.R"))

m <- 200
Time <- 200
dT <- 1
inter_T <- 1
K <- 2
sparsity <- 0.05
true_Mu <- matrix(c(2, 0.25, 0.05, 1.5), 
                  nrow = K, ncol = K, byrow = T)
# diag(true_Mu) <- c(2, 4)
# true_Mu[1, 2] <- 0.25
# true_Mu[2, 2] <- 0.75
## excitation, if used (for Hawkes)
true_B <- matrix(0, nrow = K, ncol = K, byrow = TRUE)
diag(true_B) <- 0.5
model <- "Poisson"
if(model == "Poisson") {
  true_B <- matrix(0, K, K)
}
# Pi <- c(0.2, 0.3, 0.3, 0.2)
Pi <- c(0.5, 0.5)

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
alltimes <- sampleBlockHak(Time, A, Z, Mu = true_Mu, B = true_B, lam = 1)
print("Simulated Data")

curr_n0 <- 5
m0_curr <- 50
result <- sparse_poisson(alltimes, K, n0 = curr_n0, m, m0 = m0_curr)

Mu_est <- result$est_B
## need to pass the estimated clustering also
init_tau <- matrix(0, nrow = m, ncol = K)
for(i in seq_along(result$est_clust)){
  init_tau[i, result$est_clust[i]] <- 1
}
### check the initial ARI
aricode::ARI(result$est_clust, Z)

### will need to modify to account for the decreased number
### of events also...
results_online_hom <- estimate_Poisson_init(full_data = 
                                               result$rest_events,
                                             A,
                                             m,
                                             K,
                                             Time,
                                             dT = dT,
                                             B = Mu_est,
                                             inter_T,
                                             init_tau,
                                             start = result$cut_off,
                                             is_elbo = TRUE)


plot(results_online_hom$full_ELBO, type = "l")
plot(200 * results_online_hom$wind_elbo, type = "l")
plot(results_online_hom$Cum_ELBO, type = "l")
table(Z, apply(results_online_hom$tau, 1, which.max))

## what about random init

## something not declared correctly here?

results_online <- estimate_Poisson(full_data = 
                                          as.matrix(alltimes),
                                        A,
                                        m,
                                        K,
                                        Time,
                                        dT = dT,
                                        step_size = 0.01,
                                        B = Mu_est,
                                        inter_T = 1,
                                        # init_tau,
                                        # start = 0,
                                        is_elbo = TRUE)

plot(results_online$full_ELBO, type = "l")

table(Z, apply(results_online$tau, 1, which.max))



plot(200*results_online$wind_elbo, type = "l")

## this definitely converging too fast
plot(results_online$AveELBO, type = "l") 


plot(results_online$logL, type = "l")



#### what other sort of ELBO makes sense
summary(diff(results_online_hom$full_ELBO)/max(results_online_hom$full_ELBO))


### compare this to the elbo from the corresponding batch estimator,


result_batch <- batch_estimator_hom_Poisson(
  alltimes,
  A,
  m,
  K,
  T = Time,
  itermax = 100,
  stop_eps = 0.001)

result_batch$Pi
table(Z, apply(result_batch$tau, 1, which.max))

plot(result_batch$ELBO, type = "l")

### then make a nice plot of the two of these together
elbos <- tibble(ELBO = as.numeric(result_batch$ELBO[result_batch$ELBO != 0]),
                method = "BATCH") %>% 
  bind_rows(tibble(ELBO = as.numeric(results_online_hom$full_ELBO),
                   method = "ONLINE"))

elbos %>% 
  group_by(method) %>% 
  mutate(index = row_number()) %>% 
  ggplot(aes(index, ELBO, colour = method)) +
  geom_line()


batch_data <- tibble(ELBO = as.numeric(result_batch$ELBO[result_batch$ELBO != 0]),
                     method = "BATCH") %>% 
  mutate(index = row_number(), events_seen = nrow(alltimes) * index)

init_events <- alltimes[alltimes[,3] < result$cut_off,]

online_data <- tibble(ELBO = as.numeric(results_online_hom$full_ELBO),
                      method = "ONLINE",
                      events_seen = as.numeric(results_online_hom$Cum_Events)+
                        nrow(init_events)) 

online_data


elbo_dat <- bind_rows(batch_data, online_data)

elbo_dat %>% 
  group_by(method) %>% 
  mutate(diff_elbo = ELBO - lag(ELBO)) %>% 
  ggplot(aes(events_seen, ELBO, colour = method)) +
  geom_line() +
  scale_x_log10()

elbo_dat %>% 
  ggplot(aes(events_seen, ELBO, colour = method)) +
  geom_line() 


### this sort of plot looks fine, just need to tidy it up quite a bit
### extend the batch estimate out at the end so can see it


#### thesis version plots below here.....
axis_title <- 16
axis_text_size <- 14
legend_text <- 14
line_size <- 1

## plot the full elbo first, compared to estimate reached by 
## batch 

online_data <- tibble(ELBO = as.numeric(results_online$full_ELBO),
                      method = "Online") 


batch_best_elbo <- max(result_batch$ELBO[result_batch$ELBO != 0])


batch_elbo <- tibble(index = 1:100, ELBO = batch_best_elbo,
                     method = "Batch")

plot1 <- online_data %>% 
  mutate(index = row_number()) %>% 
  bind_rows(batch_elbo) %>% 
  ggplot(aes(index, ELBO, colour = method)) +
  geom_line(size = 1) +
  geom_hline(yintercept = batch_best_elbo, colour = "red", size = 1) +
  labs(x = "Time") +
  theme(legend.title = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = axis_text_size),
        axis.title = element_text(size = axis_title),
        strip.text.x = element_text(size = axis_title),
        legend.text = element_text(size = legend_text)) 
plot1

ggsave("elbo_plot.png",
       plot1,
       path = here("output/thesis_figure/"),
       dev = "png",
       height = 5,
       width = 7)


#### number of events seen

online_data_events <- tibble(ELBO = as.numeric(results_online_hom$full_ELBO),
                             method = "Online",
                             events_seen = 
                               as.numeric(results_online_hom$Cum_Events) +
                                              nrow(init_events)) 

batch_events <- tibble(ELBO = result_batch$ELBO[result_batch$ELBO != 0],
                       method = "Batch") %>% 
  mutate(index = row_number()) %>% 
  mutate(events_seen = nrow(alltimes) * index)


##
events_data <- bind_rows(online_data_events, batch_events) %>% 
  mutate(perc_events = events_seen/nrow(alltimes))


plot2 <- events_data %>% 
  ggplot(aes(perc_events, ELBO, colour = method)) +
  geom_line(size = line_size) +
  scale_x_log10() +
  labs(x = "Events Used") +
  theme(legend.title = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = axis_text_size),
        axis.title = element_text(size = axis_title),
        strip.text.x = element_text(size = axis_title),
        legend.text = element_text(size = legend_text)) 

plot2 <- plot2 +
  scale_x_continuous(labels = scales::percent)

ggsave("elbo_plot_2.png",
       plot2,
       path = here("output/thesis_figure/"),
       dev = "png",
       height = 5,
       width = 7)

ggsave("elbo_plot_2_slides.png",
       plot2,
       path = here("output/thesis_figure/"),
       dev = "png",
       height = 3.5,
       width = 7)


### how does this compare with the normalised elbo, in terms of when/how
### it converges

plot(results_online_hom$AveELBO, type = "l")
results_online_hom$full_ELBO


tibble(full_elbo = as.numeric(results_online_hom$full_ELBO),
       avg_elbo = as.numeric(results_online_hom$AveELBO)) %>% 
  mutate(index = row_number()) %>% 
  pivot_longer(cols = full_elbo:avg_elbo) %>% 
  group_by(name) %>% 
  mutate(norm = 1 - abs(value)/max(abs(value))) %>% 
  ggplot(aes(index, norm, colour = name)) +
  geom_line()

tibble(full_elbo = as.numeric(results_online_hom$full_ELBO)) %>% 
  mutate(index = row_number()) %>% 
  mutate(norm = 2 - (full_elbo)/max((full_elbo))) %>% 
  ggplot(aes(index, norm)) +
  geom_line()


results_online_hom$full_ELBO

full_elbo <- as.numeric(results_online_hom$full_ELBO)
avg_elbo <- as.numeric(results_online_hom$AveELBO)

plot(2 - full_elbo/max(full_elbo), type = "l")

plot(1 + avg_elbo, type = "l")

norm_full_elbo <- 2 - full_elbo/max(full_elbo)

norm_avg_elbo <- 1 + avg_elbo

### then put these together

conv_elbo <- tibble(
  "Full ELBO" = norm_full_elbo,
  "Normalised ELBO" = norm_avg_elbo
) %>% 
  mutate(index = row_number()) %>% 
  pivot_longer(cols = `Full ELBO`:`Normalised ELBO`)

elbo_plot3 <- conv_elbo %>% 
  ggplot(aes(index, value, colour = name)) +
  geom_line(size = line_size) +
  labs(colour = "Metric", x = "Time", y = "ELBO") +
  theme(legend.title = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = axis_text_size),
        axis.title = element_text(size = axis_title),
        strip.text.x = element_text(size = axis_title),
        legend.text = element_text(size = legend_text)) 
  
ggsave("elbo_plot_3.png",
       elbo_plot3,
       path = here("output/thesis_figure/"),
       dev = "png",
       height = 5,
       width = 7)

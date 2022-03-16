#### analyse the results of some simulation studies
library(here)
source(here("Experiments/utils.R"))

### EXP 1, investigate convergence in time ####

exp_1_files <- list.files(path = here("Experiments/exp_results/"),
                          pattern = "exp1_time")

exp_1_results <- exp_1_files %>% 
  map(~readRDS(here("Experiments/exp_results/", .x))) %>% 
  imap(~ .x %>% map_dfr(`[`, c("clust", "time")) %>% 
         mutate(index = row_number())) %>% 
  bind_rows() %>% 
  mutate(index = row_number())
  

exp_1_results %>% 
  group_by(time) %>% 
  mutate(exact = ifelse(clust == 1, TRUE, FALSE)) %>% 
  summarise(ave_ARI = mean(clust),
            sd_ARI = sd(clust),
            prop_correct = sum(exact)/n())

exp_1_results %>% 
  group_by(time) %>% 
  mutate(time = as.factor(time)) %>% 
  ggplot(aes(time, clust)) +
  # geom_boxplot() +
  geom_point()

## similarly, look at the ELBO metrics for these simulations

# then the whole thing
elbo_dat <- exp_1_files %>% 
  map(~readRDS(here("Experiments/exp_results/", .x))) %>%
  map(~.x %>%  tibble(
    clust = map_dbl(., "clust"),
    elbo = map(., "est_elbo"),
    time = map_dbl(., "time")
  ) %>% 
    mutate(sim = row_number()) %>% 
    select(sim, clust, time, elbo) %>% 
    unnest(cols = c(elbo)) %>% 
    group_by(sim) %>% 
    mutate(wind = 1:unique(time))) %>% 
  bind_rows()
## this looks good, need to store the time window also

elbo_dat %>% 
  filter(time == 50) %>% 
  mutate(correct = ifelse(clust == 1, "Recover", "Don't")) %>% 
  ggplot(aes(wind, elbo, colour = correct)) + geom_line(aes(group = sim)) +
  labs(x = "Windows", y = "Window ELBO", colour = "Recovery") +
  geom_hline(yintercept = -0.2)

elbo_dat %>% 
  filter(time == 1000) %>% 
  filter(wind < 50) %>% 
  mutate(correct = ifelse(clust == 1, "Recover", "Don't")) %>% 
  ggplot(aes(wind, elbo, colour = correct)) + geom_line(aes(group = sim)) +
  geom_hline(yintercept = -0.05) +
  labs(x = "Windows", y = "Window ELBO", 
       title = "Comparison in initial optimisation",
       subtitle = "First 50 windows for T=500") +
  facet_wrap(~correct)

## examine the recovery of the rate matrix also ##
rate_exp_files <- list.files(path = here("Experiments/exp_results/"),
                                      pattern = "exp1_store")


rate_exp_files %>% 
  map(~readRDS(here("Experiments/exp_results/", .x))) %>% 
  flatten() %>% # to compress to a single list rather than nested
  imap(~update_list(., sim = .y)) %>% 
  map_dfr(`[`, c("sim","clust", "time", "est_B"))
  
b_data <- rate_exp_files %>% 
  map(~readRDS(here("Experiments/exp_results/", .x))) %>% 
  flatten() %>% # to compress to a single list rather than nested
  imap(~update_list(., sim = .y))


## then want this in the form of a data frame with a column for the
## matrices
## then some metric for the difference here, up to permutation
b_df <- b_data %>% {
  tibble(
    map_dfr(., magrittr::extract, c("sim", "clust", "time"))
  )
} %>% 
  mutate(sum_B = map_dbl(b_data, ~sum(.x$est_B)))
## sum entries of each and abs of difference

b_df %>% 
  mutate(diff = abs(sum_B - 3.10)/4) %>% 
  mutate(nodes = ifelse(clust == 1, "Recover", "Don't")) %>% 
  ggplot(aes(nodes, diff)) +
  geom_boxplot() +
  labs(y = "Mean sum of Difference")

b_df %>% 
  mutate(diff = abs(sum_B - 3.10)/4) %>% 
  mutate(nodes = ifelse(clust == 1, "Recover", "Don't")) %>% 
  ggplot(aes(diff)) +
  geom_histogram() +
  facet_wrap(~nodes)
  

# this looks good, but the next step of interest is to look at the
# rate of recovery of these plots over time.

all_b_files <- list.files(path = here("Experiments/exp_results/"),
                          pattern = "exp1_all")

b_ests <- all_b_files %>% 
  map(~readRDS(here("Experiments/exp_results/", .x))) %>% 
  flatten() %>% 
  imap(~update_list(., sim = .y))

sum_b <- b_ests %>% map(~ apply(.x$est_B, 3, sum) )

all_b_df <- b_ests %>% {
  tibble(
    map_dfr(., `[`, c("sim", "clust", "time")),
    sum_B = sum_b
  ) 
} %>% 
  unnest(sum_B) %>% 
  mutate(diff = abs(sum_B - 3.1)/4) %>% 
  group_by(sim) %>% 
  mutate(window = row_number()) %>% 
  ungroup()

all_b_df %>% 
  mutate(Cluster = ifelse(clust == 1, "Recover", "Don't"),
         time = as.factor(time)) %>% 
  # filter(time == 1000) %>% 
  ggplot(aes(window, diff)) + 
  geom_line(aes(group = sim, colour = Cluster)) +
  facet_wrap(~time, scales = "free_x") +
  labs(x = "Observation Period", y = "Difference in Rate",
       subtitle = "78% of Random Initialisation Recover Communities Exactly",
       title = "Recovery of Rate Matrix")

all_b_df %>% 
  mutate(ARI = ifelse(clust == 1, "Recover", "Don't")) %>% 
  group_by(ARI, time) %>% 
  distinct(sim) %>% tally()

all_b_df %>% 
  mutate(Cluster = ifelse(clust == 1, "Recover", "Don't"),
         time = as.factor(time)) %>% 
  filter(time %in% c(50, 500)) %>%
  filter(Cluster == "Recover") %>% 
  mutate(time = paste0("Time = ", time)) %>% 
  ggplot(aes(window, diff)) + 
  geom_line(aes(group = sim), colour = "Red") +
  facet_wrap(~time, scales = "free_x") +
  labs(x = "Observation Period", y = "Difference in Rate")

ggsave(filename = "exp1_rate_recov.png",
       path = here("Job_Talk_Plots"), width = 7, height = 4.5)

### EXP 2:investigate relationship with dT ####

# get files 
all_sims <- list.files(path = here("Experiments/exp_results/"),
           pattern = "exp_2")


all_output <- all_sims %>% 
  map(~readRDS(here("Experiments/exp_results/", .x)))


all_output[[1]] %>% 
  map_dfr(`[`, c("clust", "dT"))
# want this across all list elements
sims_df <- all_output %>% 
  imap(~ .x %>% map_dfr(`[`, c("clust", "dT")) %>% 
                     mutate(sim = .y)) %>% 
  bind_rows()

sims_df %>% 
  mutate(dT = as.factor(dT)) %>% 
  ggplot(aes(dT, clust)) +
  geom_boxplot() +
  labs(x = "Window Length", y = "ARI", title = "Varying window length")

sims_df %>% 
  group_by(dT) %>% 
  summarise(mean = mean(clust), sd = sd(clust),
            lower = max(mean - sd,0), upper = min(mean + sd, 1)) %>% 
  # mutate(lower = max(mean - sd,0), upper = min(mean + sd, 1)) %>% 
  ggplot(aes(dT, mean)) + 
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) +
  labs(y = "Proportion Recovered",
       subtitle = "Plus or Minus Standard Deviation",
       title = "100 Nodes, T = 100, 20 Simulations")


### EXP 3: Investigate theoretical regret rate ####
exp_3_data <- readRDS(here("Experiments/exp_results/", "exp3_200.RDS"))

exp3_tidy <- exp_3_data %>% tibble(
  sim = 1:100,
  ARI = map_dbl(., "clust"),
  regret = map(., "regret"),
  dT = list(1:200)
) %>% 
  unnest(cols = c(regret, dT)) %>% 
  select(sim, ARI, regret, dT)

exp3_tidy %>% 
  mutate(ARI = ifelse(ARI == 1, "Recover Community", "Don't Recover")) %>%
  # filter(ARI )
  filter(regret > 0) %>% 
  ggplot(aes(dT, regret)) +
  geom_line(aes(group = sim, colour = ARI), alpha = 0.4) +
  stat_function(fun = function(x) 2*sqrt(x) * log(x*1000)^2,
                colour = "black") +
  facet_wrap(~ARI, scales = "free") +
  labs(y = "Regret Rate", x = "Windows",
       title = "Comparison to theoretical rate",
       subtitle = "Ignoring the constants") +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

exp3_tidy %>% 
  group_by(sim) %>% 
  filter(regret > 0) %>% 
  group_by(ARI, sim) %>% 
  tally() %>% 
  filter(ARI == 1 & n == 200)



### EXP 4 Online Loss ####

exp_4_files <- list.files(path = here("Experiments/exp_results/"), 
                          pattern = "exp4_")

exp4_sims <- exp_4_files %>% 
  map(~readRDS(here("Experiments/exp_results/", .x))) %>% 
  flatten() %>% 
  imap(~update_list(., sim = .y)) %>% 
  map_dfr(`[`, c("online_loss", "batch_ave_loss", "sim", "clust")) 

## need to store simulation id for each, then can do a group-by and
## such

# exp4_sims %>% dim()

exp_4_tidy <- exp4_sims %>% 
  reduce(data.frame) %>% 
  as_tibble() %>% 
  rename(Batch_Loss = elt, Sim = elt.1, ARI = elt.2) %>% 
  group_by(Sim) %>% 
  mutate(Time = max(dT))


curr_time <- 200
exp_4_tidy %>% 
  filter(Time == curr_time) %>% 
  ggplot(aes(dT, Loss, colour = Z)) +
  geom_line() +
  geom_hline(aes(yintercept = Batch_Loss), 
             exp_4_tidy %>% filter(Time == curr_time)) +
  facet_wrap(~Sim, scales = "free")

## the true z not doing as well could just be an identifiability issue


## trying to verify this...
exp_4_tidy %>% 
  filter(Sim == 17)


bad_sim <- exp_4_files %>% 
  map(~readRDS(here("Experiments/exp_results/", .x))) %>% 
  flatten() %>% 
  imap(~update_list(., sim = .y)) %>% 
  map_dfr(`[`, c("z_true", "sim", "clust", "z_est")) %>% 
  filter(sim == 17) %>% 
  select(z_true, z_est)

table(bad_sim$z_true, bad_sim$z_est)

B_comp <- exp_4_files %>% 
  map(~readRDS(here("Experiments/exp_results/", .x))) %>% 
  flatten() %>% 
  imap(~update_list(., sim = .y)) %>% 
  map(`[`, c("B_ests", "sim")) 

B_comp[[17]]
# this indicates this is due to label switching alright

## pred log likelihood ###

exp4_sims_llh <- exp_4_files %>% 
  map(~readRDS(here("Experiments/exp_results/", .x))) %>% 
  flatten() %>% 
  imap(~update_list(., sim = .y)) %>% 
  map_dfr(`[`, c("sim", "pred_llh", "clust")) 

llh_tidy <- exp4_sims_llh %>% 
  reduce(data.frame) %>% 
  rename(sim = out, clust = elt) %>% 
  as_tibble()

llh_tidy %>% 
  mutate(neg_llh = - loglik) %>% 
  group_by(sim, method) %>% 
  mutate(cum_loss = cumsum(neg_llh)/dT) %>% 
  group_by(sim) %>% 
  mutate(Time = max(dT) + 1, clust = factor(clust)) %>% 
  filter(Time %in% c(50,500)) %>%
  filter(clust == 1) %>%
  # filter(sim < 10) %>%
  ggplot(aes(dT, cum_loss, colour = method)) +
  geom_line(aes(group = interaction(sim, method),
                linetype = as.factor(sim))) + 
  guides(linetype = "none") +
  # facet_wrap(vars(sim,clust), scales = "free") +
  facet_wrap(~Time, scales = "free") +
  # theme(legend.position = NULL) +
  NULL

### maybe plot the difference instead?


llh_tidy %>% 
  mutate(neg_llh = - loglik) %>% 
  group_by(sim, method) %>% 
  mutate(cum_loss = cumsum(neg_llh)/dT) %>% 
  group_by(sim) %>% 
  mutate(Time = max(dT) + 1, clust = factor(clust)) %>% 
  filter(Time %in% c(50,500)) %>%
  filter(clust == 1) %>%
  pivot_wider(id_cols = c(sim:dT, Time, clust),
              names_from = method,
              values_from = cum_loss) %>% 
  mutate(diff = online - batch) %>% 
  mutate(Time = paste0("Time = ", Time)) %>% 
  ggplot(aes(dT, diff)) +
  geom_line(aes(group = sim, colour = clust), show.legend = FALSE) + 
  # guides(linetype = "none") +
  # facet_wrap(vars(sim,clust), scales = "free") +
  facet_wrap(~Time, scales = "free") +
  labs(x = "Observation Time", y = "Difference in Loss")
 
ggsave(filename = "exp4_pred_loss_llh.png",
       path = here("Job_Talk_Plots"), width = 7, height = 4.5)

## facet these by ARI here to see

# llh_tidy %>% 
#   mutate(neg_llh = - loglik) %>% 
#   group_by(sim, method) %>% 
#   mutate(cum_loss = cumsum(neg_llh)/dT) %>% 
#   group_by(dT, method) %>% 
#   summarise(ave = mean(cum_loss), sd(cum_loss)) %>% 
#   ggplot(aes(dT, ave, colour = method)) +
#   geom_line()

### EXP 5 Changing Number of Nodes ####

exp_5_files <- list.files(path = here("Experiments/exp_results/"), 
                          pattern = "exp5_")

exp_5_data <- exp_5_files %>% 
  map(~readRDS(here("Experiments/exp_results/", .x))) %>% 
  flatten() %>% 
  imap(~update_list(., sim = .y)) %>% 
  map_dfr(`[`, c("curr_n", "clust", "sim"))


exp_5_data %>% 
  group_by(curr_n) %>% 
  summarise(mean = mean(clust), sd = sd(clust))


### EXP 6 Number of Clusters ####

exp_6_files <- list.files(path = here("Experiments/exp_results/"), 
                          pattern = "exp6_")

exp_6_data <- exp_6_files %>% 
  map(~readRDS(here("Experiments/exp_results/", .x))) %>% 
  flatten() %>% 
  imap(~update_list(., sim = .y)) %>% 
  map_dfr(`[`, c("K", "clust", "sim"))

exp_6_data %>% 
  group_by(K) %>% 
  summarise(mean = mean(clust), sd = sd(clust))


#### EXP 7 ####
exp_7_files <- list.files(path = here("Experiments/exp_results/"), 
                          pattern = "exp_7")

exp_7_data <- exp_7_files %>% 
  map(~readRDS(here("Experiments/exp_results/", .x))) %>% 
  flatten() %>% 
  imap(~update_list(., sim = .y))%>% 
  map_dfr(`[`, c("ARI", "dT", "sim", "nodes", "model"))


exp_7_data %>% 
  # group_by(nodes, model, dT) %>% 
  # summarise(mean_ari = mean(ARI), sd = sd(ARI)) %>%
  mutate(dT = as.factor(dT)) %>% 
  # filter(model == "Hawkes") %>%
  # filter(nodes == 100) %>% 
  ggplot(aes(dT, ARI)) + 
  # geom_point() +
  geom_boxplot() +
  facet_grid(cols = vars(nodes), rows = vars(model)) +
  # theme(axis.text.x = element_text(size = 4))
  scale_x_discrete(breaks = 1:5) +
  labs(x = "Window Size")




#### Online Cluster Recovery Plot ####

# which dataset is best for this? exp3?


#### Exp 9 ####

exp_9_files <- list.files(path = here("Experiments/exp_results/"), 
                          pattern = "exp9_")

exp_9_data <- exp_9_files %>% 
  map(~readRDS(here("Experiments/exp_results/", .x))) %>% 
  flatten() %>% 
  imap(~update_list(., sim = .y,
                    est_clust = apply(.x$online_clust, c(1, 3), which.max)))%>% 
  map(`[`, c("est_clust", "z_true", "sim", "clust"))



ests <- exp_9_data %>% 
  map("est_clust")

true_clust <- exp_9_data %>% 
  map("z_true")

final_est <- exp_9_data %>% 
  map("clust")

ari_ests <- map2(ests, true_clust,
                 ~ apply(.x, 2, function(x) aricode::ARI(x, .y)))
# this is what I want, now just need to convert to a tibble I can use
length(ari_ests)
names(ari_ests) <- paste0("ARI_", 1:80)


tidy_ari <- tibble(sim = 1:80, ests = ari_ests,
                   final_ari = final_est) %>% 
  rowwise() %>% 
  mutate(Time = length(ests)) %>% 
  unnest(ests) %>% 
  unnest(final_ari) %>% 
  group_by(sim) %>% 
  mutate(dT = row_number()) %>% 
  filter(Time != dT) %>% 
  mutate(Time = Time - 1) 


tidy_ari %>% 
  filter(dT %% 10 == 0) %>% 
  mutate(Time = as.factor(Time), dT = as.factor(dT)) %>% 
  # filter(final_ari == 1) %>% 
  ggplot(aes(dT, ests)) + 
  geom_boxplot() +
  facet_wrap(~Time, scales = "free")

tidy_ari %>% 
  mutate(Time = paste0("Time = ", Time)) %>% 
  mutate(Time = factor(Time, 
                       levels = paste0("Time = ",c(50, 100, 200, 500)))) %>% 
  # filter(final_ari  == 1) %>%
  ### this is cheating a bit too much
  group_by(Time, dT) %>% 
  summarise(mean_ari = mean(ests), sd = sd(ests),
            min_bar = max(mean_ari - sd, 0),
            max_bar = min(mean_ari + sd, 1)) %>% 
  filter(Time %in% c("Time = 50", "Time = 500")) %>%
  # filter(dT %% 1 == 0) %>%
  ungroup() %>% 
  ggplot(aes(dT, mean_ari)) + 
  geom_ribbon(aes(ymin = min_bar, ymax = max_bar), 
              fill = "grey70", alpha = 0.2) +
  geom_line(colour = "red") +
  facet_wrap(~Time, scales = "free") +
  labs(x = "Observation Length", y = "Estimated ARI") + ylim(0,1)

ggsave(filename = "exp1_clust_recov.png",
       path = here("Job_Talk_Plots"), width = 7, height = 4.5)

### just need to add some error bars and tidy this up a bit more


tidy_ari %>% 
  mutate(Time = paste0("Time = ", Time)) %>% 
  mutate(Time = factor(Time, 
                       levels = paste0("Time = ",
                                       c(50, 100, 200, 500)))) %>% 
  filter(dT < 25) %>% 
  ggplot(aes(dT, ests, group = sim)) +
  geom_line(colour = "red", alpha = 0.5) +
  facet_wrap(~Time, scales = "free") +
  labs(y = "ARI", x = "Observed Time") + 
  ylim(0, 1)


#### EXP 10 ####

exp_10_files <- list.files(path = here("Experiments/exp_results/"), 
                          pattern = "exp_10_")

exp_10_data <- exp_10_files %>% 
  map_dfr(~readRDS(here("Experiments/exp_results/", .x))) 

exp_10_data %>% 
  # group_by(K, model) %>% 
  # summarise(mean(ARI)) %>% 
  select(ARI:model) %>% 
  drop_na() %>%
  mutate(K = as.factor(K)) %>% 
  ggplot(aes(K, ARI)) + geom_boxplot() +
  facet_wrap(vars(model, nodes))

#### EXP 11 ####
exp_11_files <- list.files(path = here("Experiments/exp_results/"), 
                          pattern = "exp_11")

exp_11_data <- exp_11_files %>% 
  map_dfr(~readRDS(here("Experiments/exp_results/", .x))) %>% 
  mutate(sparsity = 0.25) # this is the same at the moment

exp_11_data %>% 
  group_by(nodes) %>% 
  summarise(mean(ARI), sd(ARI))

exp_11_data %>% 
  mutate(nodes = as.factor(nodes)) %>% 
  ggplot(aes(nodes, ARI)) +
  geom_boxplot()


#### Exp 12
exp_12_files <- list.files(path = here("Experiments/exp_results/"), 
                           pattern = "exp_12_")

exp_12_files %>% 
  map_dfr(~readRDS(here("Experiments/exp_results/", .x))) %>% 
  group_by(init, nodes) %>% 
  summarise(mean(ARI), sd(ARI), n())

### init is definitely making the performance worse at the moment


### Check for figure 1
fig_1_files <- list.files(path = here("Experiments/exp_results/"),
                                      pattern = "fig_1_exp_2")

fig_1_files %>% 
  map_dfr(~readRDS(here("Experiments/exp_results/", .x))) %>% 
  group_by(Method) %>% 
  summarise(mean(ARI), sd(ARI), median(ARI))
  
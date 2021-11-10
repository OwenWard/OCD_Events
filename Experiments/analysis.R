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
  summarise(mean = mean(clust), sd = sd(clust)) %>% 
  arrange(-mean)


#### EXP 3: Investigate theoretical regret rate ####
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

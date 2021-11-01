#### analyse the results of some simulation studies
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
exp_1_results <- exp_1_files %>% 
  map(~readRDS(here("Experiments/exp_results/", .x))) %>% 
  imap(~ .x %>% map_dfr(`[`, c("clust", "time", "est_elbo")) %>% 
       mutate(sim = .y)) %>% 
  bind_rows() %>% 
  mutate(index = row_number())
## not convinced this is in the format I want

a <- readRDS(here("Experiments/exp_results/", exp_1_files[3])) 

a %>% tibble(
  clust = map_dbl(., "clust"),
  elbo = map(., "est_elbo"),
  time = map_dbl(., "time")#,
  # wind = list(1:unique(time))
) %>% 
  mutate(sim = row_number()) %>% 
  select(sim, clust, time, elbo) %>% 
  unnest(cols = c(elbo)) %>% 
  group_by(sim) %>% 
  mutate(wind = 1:unique(time))

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
  labs(x = "Windows", y = "Window ELBO") +
  geom_hline(yintercept = -0.2)

elbo_dat %>% 
  filter(time == 500) %>% 
  filter(wind < 50) %>% 
  mutate(correct = ifelse(clust == 1, "Recover", "Don't")) %>% 
  ggplot(aes(wind, elbo, colour = correct)) + geom_line(aes(group = sim)) +
  geom_hline(yintercept = -0.05) +
  labs(x = "Windows", y = "Window ELBO", 
       title = "Comparison in initial optimisation",
       subtitle = "First 50 windows for T=500") +
  facet_wrap(~correct)

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
exp_3_data <- readRDS(here("Experiments/exp_results/", "exp3.RDS"))

exp3_tidy <- exp_3_data %>% tibble(
  sim = 1:100,
  ARI = map_dbl(., "clust"),
  regret = map(., "regret"),
  dT = list(1:100)
) %>% 
  unnest(cols = c(regret, dT)) %>% 
  select(sim, ARI, regret, dT)

exp3_tidy %>% 
  mutate(ARI = ifelse(ARI == 1, "Recover Community", "Don't Recover")) %>% 
  ggplot(aes(dT, regret)) +
  geom_line(aes(group = sim), alpha = 0.25) +
  stat_function(fun = function(x) sqrt(x) * log(x)^2, colour = "red") +
  facet_wrap(~ARI) +
  labs(y = "Regret Rate", x = "Windows",
       title = "Comparison to theoretical rate",
       subtitle = "Ignoring the constants") 

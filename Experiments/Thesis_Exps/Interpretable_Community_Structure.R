#### Interpretable Community Structure compared to 
#### aggregate data


### consider one day of citibike data for this?

### want an example with relatively stable community structure
### where aggregate methods don't capture something similar


library(tidyverse)
library(here)
library(lubridate)
library(kernlab)
library(osmdata)
library(ggmap)

source(here("Experiments/", "utils.R"))
source(here("functions/", "df_to_adj.R"))
source(here("functions", "pensky_fcns.R"))
source(here("functions/init_fcn.R"))

bike_data <- read_csv(here("data", "201809-citibike-tripdata.csv"))

sep_10 <- bike_data %>% 
  mutate(start_date = as_date(starttime),
         end_date = as_date(stoptime)) %>% 
  filter(start_date == "2018-09-10") %>% 
  filter(start_date == end_date)
### a monday

rm(bike_data)


sep_10 %>% 
  arrange(starttime)

start_time <- ymd_hms("2018-09-10 00:00:00", tz = "UTC")
# citi_data$minutes <- 
summary(as.numeric(sep_10$starttime - start_time))
## this is in seconds?
summary(as.numeric(sep_10$starttime - start_time)/60)
## and this is in minutes

### then tidy this up
event_data <- sep_10 %>%
  select(starttime, `start station id`, `end station id`) %>% 
  rename(start_loc = `start station id`,
         end_loc = `end station id`) %>% 
  mutate(event_time = (starttime - start_time)) %>% 
  select(-starttime) %>% 
  filter(start_loc != end_loc)

event_data %>% ggplot(aes(event_time)) +
  geom_histogram()

### then convert this to minutes or hours?

clean_data <- event_data %>% 
  mutate(minutes = as.numeric(event_time/60)) %>% 
  select(-event_time)

#### now fit spectral clustering on the counts


df <- clean_data %>% 
  rename(Send = start_loc,
         Rec = end_loc,
         Time = minutes) %>% 
  group_by(Send, Rec) %>% 
  count()

df


### get the number of nodes involved
stations <- unique(c(clean_data$start_loc, clean_data$end_loc))

sort(stations)

station_ids <- tibble(raw_id = sort(stations))

station_info <- station_ids %>% 
  mutate(new_id = row_number() - 1)

### then replace the df ids with these
num_stats <- nrow(station_info)
new_df <- df %>% left_join(station_info, by = c("Send" = "raw_id")) %>% 
  rename(new_send = new_id) %>% 
  left_join(station_info, by = c("Rec" = "raw_id")) %>% 
  rename(new_rec = "new_id") %>%
  ungroup() %>% 
  select(new_send, new_rec, n) %>% 
  rename(Send = new_send,
         Rec = new_rec) %>% 
  drop_na()

#### Clustering on count Matrix ####

A_mat <- summ_to_adj_mat(new_df, n = num_stats)

A_mat[1, ]

Lapl = diag(rowSums(A_mat)) - A_mat

eigs <- eigen(Lapl)

## then look at the eigen gap
evals <- as.numeric(eigs$values)
diffs <- diff(evals)
(diffs <- diffs[-1])
(est_K <- which.max(abs(diffs))+1)

### so 2 clusters
sc <- specc(Lapl, centers = est_K)

sc_est <- sc@.Data
### basically puts them all in one cluster, so
### no interpretation, but lets visualise anyway

sum_clust <- station_info %>% 
  mutate(clust = sc_est) %>% 
  ### then get the locations
  left_join(sep_10, by = c("raw_id" = "start station id")) %>% 
  select(raw_id:clust, `start station latitude`,
         `start station longitude`) %>% 
  rename(station_lat = `start station latitude`,
         station_long = `start station longitude`) %>% 
  distinct(raw_id, .keep_all = TRUE) %>% 
  mutate(miss_stat = is.na(station_lat))

send_nodes <- sum_clust %>% filter(miss_stat == FALSE) %>% 
  select(-miss_stat)

### some stations where no trips start in this time period

rec_nodes <- sum_clust %>% 
  filter(miss_stat == TRUE) %>% 
  left_join(sep_10, by = c("raw_id" = "end station id")) %>% 
  select(raw_id:clust, `end station latitude`,
         `end station longitude`) %>% 
  rename(station_lat = `end station latitude`,
         station_long = `end station longitude`) %>% 
  distinct(raw_id, .keep_all = TRUE)

### put back together for all stations
stats_sum_clust <- bind_rows(send_nodes,
                             rec_nodes)


stats_sum_clust %>% 
  ggplot(aes(station_long, station_lat, colour = clust)) +
  geom_point()

stats_sum_clust %>% filter(clust == 2) %>% 
  ggplot(aes(station_long, station_lat)) +
  geom_point()


bbox <- getbb("Manhattan")
# bbox[2,] <- c(40.65,40.82)

nyc_map <- get_map(bbox, maptype = "toner-background",
                   zoom = 12)
ggmap(nyc_map) +
  geom_point(data = stats_sum_clust %>% filter(clust == 1),
             aes(station_long, station_lat))


### see what proportion of trips are between these stations

## first by start station
clust_2 <- stats_sum_clust %>% 
  filter(clust == 1) %>% 
  pull(raw_id)

sep_10 %>% 
  filter(`start station id` %in% clust_2) %>% 
  group_by(`start station id`) %>% 
  tally() %>% 
  arrange(-n)

sep_10 %>% 
  group_by(`start station id`) %>% 
  tally() %>% 
  arrange(-n)

## so not the busiest in terms of origin

## then by end station
sep_10 %>% 
  filter(`end station id` %in% clust_2) %>% 
  group_by(`end station id`) %>% 
  tally() %>% 
  arrange(-n)

sep_10 %>% 
  group_by(`end station id`) %>% 
  tally() %>% 
  arrange(-n)

## so includes grand central, which has the most trips ending there,
## but no other obvious structure

### Adaptive Clustering of Pensky and Zhang

## 24 windows (for each hour)

## need to change this to new id's
station_info

event_data <- clean_data %>% 
  left_join(station_info, by = c("start_loc" = "raw_id")) %>% 
  rename(V1 = new_id) %>% 
  left_join(station_info, by = c("end_loc" = "raw_id")) %>% 
  rename(V2 = new_id, V3 = minutes) %>%
  select(V1, V2, V3)


A_mats <- event_to_mat_seq(event_data,
                           Total_time = 60*24,
                           window_size = 60,
                           n = num_stats)

dim(A_mats) # looks good

A_mats[A_mats > 0] <- 1

l <- 4
m_pz <- 4
m0 <- 1
l0 <- 2

A_pz <- pz_estimator_3(A = A_mats, time = 60*24,
                       l0 = l0, m0 = m0,
                       m = m_pz)

### then look at eigen gap of this
A_pz

eigs_pz <- eigen(A_pz)
evals_pz <- as.numeric(eigs_pz$values)
diffs_pz <- diff(evals_pz)
diffs_pz <- diffs_pz[-1]
(est_K_pz <- which.max(abs(diffs_pz)) + 1)

## this suggests K = 2 also

clust_pz <- specc(A_pz, centers = est_K_pz)

table(clust_pz@.Data, sc_est)
## so no overlap in the small cluster, which is fine

## visualise these clusters

pz_clust_send <- station_info %>% 
  mutate(clust = clust_pz) %>% 
  ### then get the locations
  left_join(sep_10, by = c("raw_id" = "start station id")) %>% 
  select(raw_id:clust, `start station latitude`,
         `start station longitude`) %>% 
  rename(station_lat = `start station latitude`,
         station_long = `start station longitude`) %>% 
  distinct(raw_id, .keep_all = TRUE) %>% 
  mutate(miss_stat = is.na(station_lat))

pz_send_nodes <- pz_clust_send %>%
  filter(miss_stat == FALSE) %>% 
  select(-miss_stat)

### some stations where no trips start in this time period

pz_rec_nodes <- pz_clust_send %>% 
  filter(miss_stat == TRUE) %>% 
  left_join(sep_10, by = c("raw_id" = "end station id")) %>% 
  select(raw_id:clust, `end station latitude`,
         `end station longitude`) %>% 
  rename(station_lat = `end station latitude`,
         station_long = `end station longitude`) %>% 
  distinct(raw_id, .keep_all = TRUE) 

### put back together for all stations
pz_clust <- bind_rows(pz_send_nodes,
                      pz_rec_nodes) 

table(pz_clust$clust)

pz_clust %>% 
  filter(clust == 2) %>% 
  ggplot(aes(station_long, station_lat)) +
  geom_point()

bbox <- getbb("Manhattan")
bbox[2,] <- c(40.65,40.82)

nyc_map_2 <- get_map(bbox, maptype = "toner-background",
                   zoom = 12)

ggmap(nyc_map_2) +
  geom_point(data = pz_clust,
             aes(station_long, station_lat, colour = as.factor(clust))) +
  labs(x = element_blank(), y = element_blank()) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())


### see where these are
locs <- pz_clust %>% 
  filter(clust == 2) %>% 
  select(station_lat, station_long)

ind <- 1
cat(as.numeric(locs[ind, 1]), as.numeric(locs[ind, 2]))
ind <- 2


### so these don't make any sort of sense

## confirm that the counts aren't crazy
clust_1_pz <- pz_clust %>% 
  filter(clust == 1) %>% 
  pull(raw_id)

clust_1_pz

sep_10 %>% 
  filter(`start station id` %in% clust_1_pz) %>% 
  group_by(`start station id`) %>% 
  tally()

sep_10 %>% 
  group_by(`start station id`) %>% 
  tally() %>% 
  arrange(-n)

sep_10 %>% 
  filter(`start station id` %in% clust_1_pz) %>% 
  group_by(`end station id`) %>% 
  tally() %>% 
  arrange(-n)

sep_10 %>% 
  group_by(`end station id`) %>% 
  tally() %>% 
  arrange(-n)

### trips ending at these stations

sep_10 %>% 
  filter(`end station id` %in% clust_1_pz) %>% 
  group_by(`end station id`) %>% 
  tally()

sep_10 %>% 
  group_by(`end station id`) %>% 
  tally() %>% 
  arrange(-n) %>% 
  print(n = 15)

sep_10 %>% 
  filter(`end station id` %in% clust_1_pz) %>% 
  group_by(`start station id`) %>% 
  tally() %>% 
  arrange(-n)

### so the PZ clusters are consistent across time, even 
### if the ones from the count matrix aren't

### Fit OCD Method to this Data #####

K <- 3
result <- sparse_poisson(alltimes = event_data,
                         K,
                         n0 = 180,
                         m0 = 100,
                         m = num_stats)

### need to construct A
event_data <- event_data %>% drop_na()

A = list()

for(i in 1:nrow(event_data)){
  j = event_data$V1[i]
  k = event_data$V2[i]
  #print(j)
  A[[j+1]] = j
  A[[k+1]] = k
  #A[[j+1]] = c(A[[j+1]],emails$Rec[i])
}

for(i in 1:nrow(event_data)){
  j = event_data$V1[i]
  k = event_data$V2[i]
  #print(j)
  #A[[j+1]] = c(A[[j+1]],j)
  #A[[k+1]] = c(A[[k+1]],k)
  A[[j+1]] = c(A[[j+1]], event_data$V2[i])
}


A_test = lapply(A,unique)


## Fitting Inhomogeneous Poisson
dT <- 15
H <- 2
MuA <- array(runif(K * K * H), c(K, K, H))
window <- 30
### construct tau
init_tau <- matrix(0, nrow = num_stats, ncol = K)
for(i in seq_along(result$est_clust)){
  init_tau[i, result$est_clust[i]] <- 1
}

Time <- 60*24

results_online_inpois <- nonhomoPois_est_init(alltimes = 
                                                as.matrix(result$rest_events),
                                              A = A_test,
                                              m = num_stats,
                                              K,
                                              H,
                                              window,
                                              Time,
                                              dT,
                                              gravity = 0.001,
                                              MuA, 
                                              init_tau,
                                              start = result$cut_off)

## so this runs here
## does it give a stable result?
clust <- apply(results_online_inpois$tau, 1, which.max)
table(clust)

## this looks half reasonable
results_online_inpois$MuA


### what about hom Poisson?
Mu_est <- result$est_B
## need to pass the estimated clustering also
init_tau <- matrix(0, nrow = num_stats, ncol = K)
for(i in seq_along(result$est_clust)){
  init_tau[i, result$est_clust[i]] <- 1
}
results_online_init <- estimate_Poisson_init(full_data = 
                                               as.matrix(result$rest_events),
                                             A = A_test,
                                             m = num_stats,
                                             K,
                                             Time,
                                             dT = dT,
                                             B = Mu_est,
                                             inter_T = 1,
                                             init_tau,
                                             start = result$cut_off,
                                             is_elbo = FALSE)

results_online_init$Pi

table(apply(results_online_init$tau, 1, which.max))





##### Remove night time Data ####

day_data <- event_data %>% 
  mutate(V3 = V3 - 360) %>% 
  filter(V3 > 0) %>% 
  filter(V3 < 900)


### so now total = 60*18 = 1080
K <- 3
result <- sparse_poisson(alltimes = day_data,
                         K,
                         n0 = 60,
                         m0 = 100,
                         m = num_stats)

A = list()

for(i in 1:nrow(day_data)){
  j = day_data$V1[i]
  k = day_data$V2[i]
  #print(j)
  A[[j+1]] = j
  A[[k+1]] = k
  #A[[j+1]] = c(A[[j+1]],emails$Rec[i])
}

for(i in 1:nrow(day_data)){
  j = day_data$V1[i]
  k = day_data$V2[i]
  #print(j)
  #A[[j+1]] = c(A[[j+1]],j)
  #A[[k+1]] = c(A[[k+1]],k)
  A[[j+1]] = c(A[[j+1]], day_data$V2[i])
}


A_test = lapply(A, unique)

dT <- 15
H <- 2
MuA <- array(runif(K * K * H), c(K, K, H))
window <- 15
### construct tau
init_tau <- matrix(0, nrow = num_stats, ncol = K)
for(i in seq_along(result$est_clust)){
  init_tau[i, result$est_clust[i]] <- 1
}

Time <- 60*24

results_online_inpois <- nonhomoPois_est_init(alltimes = 
                                                as.matrix(result$rest_events),
                                              A = A_test,
                                              m = num_stats,
                                              K,
                                              H,
                                              window,
                                              Time,
                                              dT,
                                              gravity = 0.001,
                                              MuA, 
                                              init_tau,
                                              start = result$cut_off)

table(apply(results_online_inpois$tau, 1, which.max))

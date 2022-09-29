#### debugging the init function in the dense Poisson case
### trying to identify the reason it doesn't work

library(here)
library(kernlab)

source(here("Experiments/", "utils.R"))
source(here("functions/init_fcn.R"))
source(here("functions/", "df_to_adj.R"))


### simulate some dense Poisson Data

Time <- 200
dT <- 1
inter_T <- 1
K <- 2
m_vec <- c(100, 200, 400)
sparsity <- 0.5 # prop of edges which can have events


model <- "Poisson"


results <- list()
m <- 100

# for(exp_num in seq_along(m_vec)) {
#   dT <- 1
curr_dt_sims <- tibble()
cat("Current K:", K, "\n")
cat("Current m:", m, "\n")
cat(model, "\n")
## baseline rate of the process
true_Mu <- matrix(0.05, 
                  nrow = K, ncol = K, byrow = T)
diag(true_Mu) <- c(rep(0.5, K-1), 1)
## excitation, if used (for Hawkes)
true_B <- matrix(0, nrow = K, ncol = K, byrow = TRUE)
diag(true_B) <- 0.5
if(model == "Poisson") {
  true_B <- matrix(0, K, K)
}
Pi <- matrix(1/K, 1, K)

Z <- sample(0:(K-1), size = m, prob = Pi, replace = TRUE)
# then generate A
A <- list()
for(i in 1:m){
  # could sample these here with SBM structure...
  num_edge = m * sparsity
  edge <- sample(m, num_edge) - 1
  ### remove self edges in here instead
  edge <- sort(edge[edge!= i-1])
  
  # edge <- sort(edge)
  A[[i]] <- edge
}
alltimes <- sampleBlockHak(Time, A, Z, Mu = true_Mu, B = true_B, lam = 1)


result <- dense_poisson(alltimes, K, n0 = 30)
aricode::ARI(result$est_clust, Z)
### get the cut off from this

cutoff <- result$cut_off

### check the performance of just using spectral clustering on the events
### before the cut off


init_events <- alltimes[alltimes[,3] < cutoff, ]

df <- as_tibble(init_events) %>% 
  rename(Send = V1,
         Rec = V2,
         Time = V3) %>% 
  filter(Send != Rec) %>% 
  group_by(Send, Rec) %>% 
  count()

A_mat <- summ_to_adj_mat(df, n = m)

dim(A_mat)

sc <- specc(A_mat, centers = 2)
sc_est <- sc@.Data
(sc_ari <- aricode::ARI(Z, sc_est))

## so spectral clustering recovers them exactly

(init_ari <- aricode::ARI(result$est_clust, Z))

## mean while our init procedure doesn't recover them at all

### what about if we use all the data in the init procedure

result <- dense_poisson(alltimes, K, n0 = 50)

(init_ari <- aricode::ARI(result$est_clust, Z))

### still doesn't work as well


### debugging the components of this init function directly


tidy_events <- as_tibble(alltimes) %>% 
  rename(send = V1, rec = V2, time = V3)
N <- nrow(tidy_events)
C <- 1
n0 <- C * log(N)
### set it for now
n0 <- 200


init_events <- tidy_events %>% 
  filter(time <= n0) %>% 
  filter(send != rec)
remaining_events <- alltimes[alltimes[,3]> n0, ]

out_events <- init_events %>% 
  group_by(send) %>% 
  count() %>% 
  ungroup() %>% 
  rename(out = n)

max_degree <- out_events %>% 
  slice_max(out, n = 1) %>% 
  pull(send)

est_int <- init_events %>% 
  filter(send == max_degree) %>% 
  group_by(rec) %>% 
  count() %>% 
  mutate(est_lam = n/n0) %>% 
  select(rec, est_lam)

est <- kmeans(est_int$est_lam, centers = K)
clust_ests <- tibble(node = est_int$rec, clust = est$cluster)

neighbours <- init_events %>% 
  filter(send == max_degree) %>% 
  ### check this giving a warning too
  distinct(rec) %>% 
  arrange(rec) %>% 
  pull(rec)


neighbour_events <- init_events %>% 
  filter(send %in% neighbours) %>% 
  filter(rec %in% neighbours) %>% 
  group_by(send, rec) %>% 
  count() %>% 
  ## then join the cluster assignments here
  left_join(clust_ests, by = c("send" = "node")) %>% 
  rename(send_clust = clust) %>% 
  left_join(clust_ests, by = c("rec" = "node")) %>% 
  rename(rec_clust = clust)

init_B <- matrix(NA, nrow = K, ncol = K)

for(k1 in 1:K){
  for(k2 in 1:K){
    curr_est <- neighbour_events %>% 
      filter(send_clust == k1) %>% 
      filter(rec_clust == k2) %>% 
      ungroup() %>% 
      slice_sample(n = 1) %>% 
      mutate(rate = n/n0) %>% 
      pull(rate)
    ### check if empty here
    if(identical(curr_est, numeric(0))){
      init_B[k1, k2] <- 0
    }
    else{
      init_B[k1, k2] <- curr_est 
    }
  }
}

### this looks like a reasonable estimate
init_B
true_Mu

init_group <- rep(NA, m)

dist_matrix <- matrix(NA, nrow = m, ncol = K)
b_ests <- matrix(NA, nrow = m, ncol = K)

# init_B <- true_Mu
### just to check

for(i in (1:m)){
  ## 
  # print(i)
  curr_neigh <- init_events %>% 
    filter(send == i-1) %>% 
    group_by(send, rec) %>% 
    summarise(count = n()) %>% 
    mutate(rate = count/n0) %>% 
    pull(rate)
  # print(length(curr_neigh))
  # then k means on these estimates
  ## check number of neighbours here too
  curr_est <- kmeans(curr_neigh, centers = K)
  curr_center <- curr_est$centers
  # print("K Means Center")
  # print(curr_center)
  
  # print(init_B)
  ## then see which row of init_B this is closest to
  
  ### this to check if there is an issue with permutations
  # perm_B <- rbind(init_B, init_B[,c(2,1)])
  # locs <- c(1, 2, 1, 2)
  
  ### this a way to try get around it?
  # curr_center <- sort(curr_center)
  perm_B <- init_B
  locs <- order(curr_center)
  
  b_ests[i,] <- as.vector(curr_center)
  ### ---- 
  dists <- apply(perm_B, 1, function(x) 
    ### change from rows to columns maybe?
    dist(rbind(x, as.vector(curr_center))))
  dist_matrix[i,] <- dists
  # print("Distances")
  # print(dists)
  # print(which.min(dists))
  # init_group[i] <- locs[which.min(dists)]
  # print("=======")
}


### this appears to be poor, 
### still doesn't work
table(init_group, Z)
aricode::ARI(init_group, Z)


clust_info <- tibble(nodes = 0:(m-1), clust = init_group)

clust_events <- init_events %>% 
  left_join(clust_info, by = c("send" = "nodes")) %>% 
  rename(clust_send = clust) %>% 
  left_join(clust_info, by = c("rec" = "nodes")) %>% 
  rename(clust_rec = clust)

updated_B <- init_B

for(k1 in 1:K){
  for(k2 in 1:K){
    new_est <- clust_events %>% 
      filter(clust_send == k1) %>% 
      filter(clust_rec == k2) %>% 
      group_by(send, rec) %>% 
      summarise(num_events = n()) %>% 
      ungroup() %>% 
      ### fitting a common Poisson to each of these pairs
      summarise(est_rate = sum(num_events)/ (n()*n0)) %>% 
      pull(est_rate)
    
    if(identical(new_est, numeric(0))){
      updated_B[k1, k2] <- 0
    }
    else{
      updated_B[k1, k2] <- new_est 
    }
  }
}

updated_B
### this is much further away than 
init_B



#### update
## what about passing the estimates from the init_B 
## into this somehow

kmeans(curr_neigh, centers = c(0.5, 1))



curr_center <- c(0.072, 0.97454)
init_B

apply(init_B, 1, function(x) 
  dist(rbind(x, as.vector(curr_center))))

dist(rbind(init_B[,1], curr_center))
dist(rbind(init_B[,2], curr_center))

#### Enron Data in R ###
library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(lubridate)
sourceCpp("C:/Users/owenw/Desktop/Online_Point_Process/onlineblock.cpp")
#set.seed(200)
library(igraph)
library(igraphdata)

data("enron")
enron[[3]]

enron[[9]][[4]]$LDC_topic

enron_node_info <- vertex.attributes(enron)

names <- enron_node_info$Name

enron_attrib <- get.edge.attribute(enron)
enron_attrib$Time
enron_edges <- ends(enron,E(enron))
enron_edges

enron_data <- cbind(enron_edges,enron_attrib$Time)

dim(enron_data)


enron_df <- tibble(send=enron_data[,1],
                        rec = enron_data[,2],
                        time = enron_data[,3])
enron_df %>% head()

lower_time <- "2000-04-27 00:00:00 UTC"
lower_time <- ymd_hms(lower_time)

enron_df <- enron_df %>% mutate(send = as.numeric(send),
                  rec = as.numeric(rec),
                  time = ymd_hms(time),
                  year = year(time)) %>%
  filter(time > lower_time) %>%
  filter(send != rec) %>%
  distinct(,.keep_all = TRUE)

na_locs <- which(names =="NA")

enron_df %>% head()

# select the data for the algorithm

# min_date <- min(enron_df$time)
# 
# enron_dat <- enron_df %>%
#   mutate(dur = time-min_date,
#          send_id = send -1,
#          rec_id = rec -1) %>%
#   select(send_id,rec_id,dur)
# 
# enron_dat %>% head()

# remove the people who are unknown
names_keep <- names[-na_locs]
enron_node_info$Note[-na_locs]

id_keep <- 1:length(names)
id_keep <- id_keep[-na_locs]

# recode based on this
enron_dat <- enron_df %>%
  filter(send %in% id_keep & rec %in% id_keep) 

enron_dat %>% head()



enron_dat <- enron_dat %>% 
  mutate(days = as.numeric(dur/(24*3600)))
# this is in days now

summary(enron_dat$days)


### then construct A
A = list()

for(i in 1:nrow(enron_dat)){
  j = enron_dat$send_id[i]
  k = enron_dat$rec_id[i]
  #print(j)
  A[[j+1]] = j
  A[[k+1]] = k
  #A[[j+1]] = c(A[[j+1]],emails$Rec[i])
}

#init_A <- A

for(i in 1:nrow(enron_dat)){
  j = enron_dat$send_id[i]
  k = enron_dat$rec_id[i]
  #print(j)
  #A[[j+1]] = c(A[[j+1]],j)
  #A[[k+1]] = c(A[[k+1]],k)
  A[[j+1]] = c(A[[j+1]],enron_dat$rec_id[i])
}


m <- length(unique(c(enron_dat$send_id,enron_dat$rec_id)))

A_test = lapply(A,unique)

# 72 and 118 are empty here
# but these will be messed up now afterwards
A_test[[72]] <- NULL
A_test[[118]] <- NULL


rm(A)


K <- 5
#m <- 184
Pi = rep(1/K,K)
B = matrix(runif(K*K),K,K)
#diag(B) = rnorm(K,mean = 0.5, sd = 0.1)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)

dT <- 20 # this is for a single hour now
inter_T <- 10


Time = 1320 # need to change this if the data changes

enron_final <- enron_dat[,c(1,2,4)]



online_pois <- 
  estimate_Poisson(full_data = as.matrix(enron_final),
                            A_test,m,K,
                            T=Time,
                            dT,B,
                            tau,Pi,S,inter_T,is_elbo = TRUE)
online_pois$Pi

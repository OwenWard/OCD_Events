library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(lubridate)
library(fcd)

sourceCpp("C:/Users/owenw/Desktop/Online_Point_Process/onlineblock.cpp")

london_data <- read_csv("C:/Users/owenw/Google Drive/Tian/Oral/data/London_166JourneyDataExtract12Jun2019-18Jun2019.csv")


london_data %>% head()


london_bike <- london_data %>% dplyr::select(`StartStation Id`,`EndStation Id`,
                                      `Start Date`,
                                      `End Date`)


summary(london_bike$`EndStation Id`)

start_time <- ymd_hms("2019-06-12 00:00:00")
london_bike_dat <- london_bike %>% 
  mutate(Time = as.POSIXct(`Start Date`,
                           format = "%d/%m/%Y %H:%M") - start_time ) %>%
  arrange(Time) %>%
  mutate(start_stat = as.numeric(as.factor(`StartStation Id`))-1,
         end_stat = as.numeric(as.factor(`EndStation Id`))-1) %>%
  filter(start_stat != end_stat)

ggplot(london_bike_dat) + geom_histogram(aes(Time))

# then get m for this
m <- length(unique(c(london_bike_dat$start_stat,london_bike_dat$end_stat)))
m


model_dat <- london_bike_dat %>% dplyr::select(start_stat,end_stat,Time)

A = list()

for(i in 1:nrow(model_dat)){
  j = model_dat$start_stat[i]
  k = model_dat$end_stat[i]
  #print(j)
  A[[j+1]] = j
  A[[k+1]] = k
  #A[[j+1]] = c(A[[j+1]],emails$Rec[i])
}

for(i in 1:nrow(model_dat)){
  j = model_dat$start_stat[i]
  k = model_dat$end_stat[i]
  #print(j)
  #A[[j+1]] = c(A[[j+1]],j)
  #A[[k+1]] = c(A[[k+1]],k)
  A[[j+1]] = c(A[[j+1]],model_dat$end_stat[i])
}


A_test = lapply(A,unique)
rm(A)


K <- 5
#m <- 184
Pi = rep(1/K,K)
B = matrix(runif(K*K),K,K)
#diag(B) = rnorm(K,mean = 0.5, sd = 0.1)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)

dT <- 1 # this is for a single week now
inter_T <- 1


Time = 7*24 # need to change this if the data changes

model_dat$start_stat = as.numeric(model_dat$start_stat)
model_dat$end_stat = as.numeric(model_dat$end_stat)
#model_dat$Time = as.numeric(model_dat$Time)/24
model_dat$Time = as.numeric(model_dat$Time) # in hours instead of days
online_pois <- estimate_Poisson(full_data = as.matrix(model_dat),
                                #full_data = as.matrix(sample_data),
                                A_test,m,K,
                                T=Time,
                                dT,B,
                                tau,Pi,S,inter_T,is_elbo = TRUE)

online_pois$Pi
online_pois$AveELBO[length(online_pois$AveELBO)]*nrow(model_dat)
online_pois$logL[length(online_pois$logL)]*nrow(model_dat)

batch_pois <-batch_estimator_hom_Poisson(as.matrix(model_dat),
                                         A_test,m,K,
                                         T=Time,
                                         B,tau,               
                                         itermax = 50,
                                         stop_eps = 0.01)

tail(batch_pois$ELBO,1)
batch_pois$Pi
#tail(batch_pois$LL,1)



online_pois_init <- estimate_Poisson(full_data = as.matrix(model_dat),
                                     #full_data = as.matrix(sample_data),
                                     A_test,m,K,
                                     T=Time,
                                     dT,B=batch_pois$B,
                                     tau=batch_pois$tau,
                                     Pi=batch_pois$Pi,S,inter_T,is_elbo = TRUE)

online_pois_init$Pi
tail(online_pois_init$logL*nrow(model_dat),1)
tail(online_pois$logL*nrow(model_dat),1)
tail(batch_pois$ELBO,1)
tail(online_pois_init$AveELBO*nrow(model_dat),1)

### spectral clustering

counts <- model_dat %>% group_by(start_stat,end_stat) %>% tally()

A <- matrix(0,nrow = m,ncol=m)
for(i in 1:nrow(counts)){
  row <- counts$start_stat[i]+1
  col <- counts$end_stat[i] + 1
  A[row,col] <- counts$n[i] 
}


A_norm <- laplacian(A,normalised = T)
#eigs <- eigen(A)
#plot(sort(real(eigs$values)),type = "l")


eigs_norm <- eigen(A_norm,symmetric = TRUE)
plot(eigs_norm$values[1:10])
# this seems to indicate one cluster

bin_A <- A
bin_A[bin_A>0] = 1

norm_Bin <- laplacian(bin_A,normalised = T)
eigs_bin <- eigen(norm_Bin,symmetric = TRUE)
plot(eigs_bin$values[1:10])
# this doesn't give anything



#### non homo poisson
window = 1
K <- 3 
H <- 7
dT = 1 # 
MuA_start = array(0,c(K,K,H))
tau_start = matrix(runif(m*K),m,K)
B_start = matrix(0,K,K)


system.time(batch <- 
              batch_nonhomoPois_estimator(as.matrix(model_dat),
                                          A_test,m,K,H,window,
                                          T = Time,dT,
                                          gravity = 0.001,MuA_start,tau_start,
                                          itermax = 100,stop_eps = 0.01 ))
batch$Pi


system.time(non_hom_pois_est <- 
              nonhomoPois_estimator(as.matrix(model_dat),
                                    A_test,m,K,H,window,
                                    T = Time,dT,
                                    gravity = 0.001,
                                    MuA_start,# = batch$MuA,
                                    tau_start))# = batch$tau) )

non_hom_pois_est$Pi


batch_hawkes_nh <- batch_nonhomoHak_estimator(as.matrix(model_dat),
                                              A_test,m,K,H,window,
                                              T = Time,dT,
                                              lam = 1.5,
                                              gravity = 0.001,
                                              B_start,
                                              MuA_start,tau_start,
                                              itermax = 100,stop_eps = 0.01 )

online_hawkes_nh <- nonhomoHak_estimator_eff_revised(as.matrix(model_dat),
                                                     A_test,m,K,H,window,
                                                     T = Time,dT=0.25,
                                                     lam = 1.5,
                                                     gravity = 0.001,
                                                     B_start = batch_hawkes_nh$B,
                                                     MuA_start = batch_hawkes_nh$MuA,
                                                     tau_start = batch_hawkes_nh$tau)
online_hawkes_nh$Pi

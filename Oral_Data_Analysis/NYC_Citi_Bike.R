library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(lubridate)
sourceCpp("C:/Users/owenw/Desktop/Online_Point_Process/onlineblock.cpp")
#set.seed(200)



nyc_data <- read_csv("C:/Users/owenw/Google Drive/Tian/Oral/data/201907-citibike-tripdata.csv.zip")

# can use starttime or stoptime here
citi_data <- nyc_data %>% dplyr::select(`start station id`,`end station id`,stoptime,`start station latitude`,
                                 `start station longitude`,
                                 `end station latitude`,
                                 `end station longitude`) %>% arrange(stoptime)

citi_data %>% ggplot(aes(stoptime)) + geom_histogram()

rm(nyc_data)

summary(citi_data$stoptime)
# so data over approximately one month
# want to convert this into days
start_time <- ymd_hms("2019-07-01 00:00:00", tz = "UTC")
citi_data$minutes <- citi_data$stoptime-start_time

head(citi_data$minutes)



# this will actually be hours now
citi_data <- citi_data %>% mutate(days = as.numeric(minutes/60)) %>%
  filter(days < 7*24)



# sometimes this minutes is in minutes and sometimes this
# is in seconds.... 
# not sure how you can specify it. check with
dim(citi_data)


citi_data %>% group_by(`start station id`,`end station id`) %>% 
  tally() %>% arrange(desc(n)) %>% 
  filter(`start station id`!=`end station id`) 

#citi_data %>% ggplot(aes(days)) + geom_histogram()
#citi_data %>% ggplot(aes(stoptime)) + geom_histogram()

# so a pretty small proportion of these trips are between the same
# stations

citi_data <- citi_data  %>%
  filter(`start station id`!= `end station id`) %>%
  filter(`end station id`!= 3791) # jersey city

bike_data <- citi_data %>% drop_na() %>% 
  dplyr::select(`start station id`,`end station id`,days)
# then convert the stations to be numeric factors which start at 0


bike_data$`start station id` <- as.numeric(as.factor(bike_data$`start station id`))-1
bike_data$`end station id` <- as.numeric(as.factor(bike_data$`end station id`))-1

bike_data %>% group_by(`start station id`,`end station id`) %>% tally()

m <- length(unique(c(bike_data$`start station id`,
                     bike_data$`end station id`)))

citi_data$start_stat <- bike_data$`start station id` +1
citi_data$end_stat <- bike_data$`end station id`+1
# can use this to join then after





A = list()

for(i in 1:nrow(bike_data)){
  j = bike_data$`start station id`[i]
  k = bike_data$`end station id`[i]
  #print(j)
  A[[j+1]] = j
  A[[k+1]] = k
  #A[[j+1]] = c(A[[j+1]],emails$Rec[i])
}

for(i in 1:nrow(bike_data)){
  j = bike_data$`start station id`[i]
  k = bike_data$`end station id`[i]
  #print(j)
  #A[[j+1]] = c(A[[j+1]],j)
  #A[[k+1]] = c(A[[k+1]],k)
  A[[j+1]] = c(A[[j+1]],bike_data$`start station id`[i])
}


A_test = lapply(A,unique)
rm(A)





K <- 3
#m <- 184
Pi = rep(1/K,K)
B = matrix(runif(K*K),K,K)
#diag(B) = rnorm(K,mean = 0.5, sd = 0.1)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)

dT <- 1 # this is for a single day now
inter_T <- 1


Time = 7*24 # need to change this if the data changes

online_pois <- estimate_Poisson(full_data = as.matrix(bike_data),
                                #full_data = as.matrix(sample_data),
                                A_test,m,K,
                                T=Time,
                                dT,B,
                                tau,Pi,S,inter_T,is_elbo = TRUE)
online_pois$Pi
online_pois$AveELBO[length(online_pois$AveELBO)]*dim(bike_data)[1]
online_pois$logL[length(online_pois$logL)]*dim(bike_data)[1]


plot(online_pois$AveELBO,type="l")



start_stations <- citi_data %>% 
  dplyr::select(`start station latitude`,`start station longitude`,
         start_stat) %>% 
  distinct(start_stat,.keep_all = T) %>%
  dplyr::select(stat_id = start_stat,stat_lat = `start station latitude`,
         stat_lon = `start station longitude`)

end_stations <- citi_data %>% 
  dplyr::select(`end station latitude`,`end station longitude`,
         end_stat) %>% 
  distinct(end_stat,.keep_all = T) %>%
  dplyr::select(stat_id = end_stat,stat_lat = `end station latitude`,
         stat_lon = `end station longitude`)

stations <- bind_rows(start_stations,end_stations) %>%
  distinct(stat_id,.keep_all = T)


est_Z_online <- apply(online_pois$tau,1,which.max)
station_clusts <- tibble(station = 1:m,clust = est_Z_online)

nyc_map_data <- stations %>% left_join(station_clusts, 
                                       by = c("stat_id" = "station"))

stations %>% left_join(station_clusts, 
                       by = c("stat_id" = "station")) %>%
  ggplot(aes(stat_lon,
             stat_lat,color=as.factor(clust))) + 
  geom_point()



library(osmdata)
#library(sf)
library(ggmap)

bbox <- getbb("Manhattan")
bbox[2,] <- c(40.65,40.82)

nyc_map <- get_map(bbox, maptype = "toner-background",
                   zoom = 12)

ggmap(nyc_map) + geom_point(data = nyc_map_data,
                            aes(stat_lon,                                                    stat_lat,color=as.factor(clust))) +
  theme_nothing()



clust_data <- nyc_map_data %>% filter(clust == 3)
ggmap(nyc_map) + geom_point(data = clust_data,
                            aes(stat_lon,                                                    stat_lat)) +
  theme_nothing() 


# this useful if one interesting small cluster
library(revgeo)
revgeo(longitude = clust_data$stat_lon,
       latitude = clust_data$stat_lat,
       provider = 'photon',output = 'frame') 


### batch poisson #####

batch_pois <-batch_estimator_hom_Poisson(as.matrix(bike_data),
                                         A_test,m,K,
                                         T=Time,
                                         B,tau,               
                                         itermax = 50,
                                         stop_eps = 0.01)
batch_pois$Pi
plot(batch_pois$ELBO,type="l")

est_Z_batch <- apply(batch_pois$tau,1,which.max)



online_init <- estimate_Poisson(as.matrix(bike_data),
                                A_test,m,K,
                                T = Time, dT,
                                B = batch_pois$B,
                                tau = batch_pois$tau,Pi,S,inter_T,is_elbo = T)
online_init$Pi


online_pois$AveELBO[length(online_pois$AveELBO)]*dim(bike_data)[1]
batch_pois$ELBO[length(batch_pois$ELBO)]

online_pois$logL[length(online_pois$logL)]*dim(bike_data)[1]
batch_pois$LL[length(batch_pois$LL)]

station_clusts <- tibble(station = 1:m,clust = est_Z_batch)

stations %>% left_join(station_clusts, 
                       by = c("stat_id" = "station")) %>%
  ggplot(aes(stat_lon,
             stat_lat,color=as.factor(clust))) + 
  geom_point()

nyc_map_data <- stations %>% left_join(station_clusts, 
                                       by = c("stat_id" = "station"))


ggmap(nyc_map) + geom_point(data = nyc_map_data,
                            aes(stat_lon,                                                    stat_lat,color=as.factor(clust))) +
  theme_nothing()



clust_data <- nyc_map_data %>% filter(clust == 2)
revgeo(longitude = clust_data$stat_lon,
       latitude = clust_data$stat_lat,
       provider = 'photon',output = 'frame') 



#### Hawkes Process ####

K <- 3
Pi = rep(1/K,K)
Mu = matrix(runif(K*K),K,K)
B = matrix(runif(K*K),K,K)
#diag(B) = rnorm(3,mean = 1, sd = 0.1)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)


results_hawkes <-
  online_estimator_eff_revised(as.matrix(bike_data),
                               A_test, 
                               m, K,
                               Time,
                               dT, 
                               lam = 1.75,
                               B, Mu,
                               tau,
                               inter_T = 1, is_elbo = T)

results_hawkes$Mu
results_hawkes$B
est_Z_hawkes <- apply(results_hawkes$tau,1,which.max)
results_hawkes$Pi
results_hawkes$elbo[length(results_hawkes$elbo)]*dim(bike_data)[1]


batch_hawkes <- batch_estimator(as.matrix(bike_data),
                                A_test, 
                                m, K,
                                Time,
                                dT, 
                                lam = 1,
                                B, Mu,
                                tau,
                                itermax = 50,
                                stop_eps = 0.01)

batch_hawkes$Mu
batch_hawkes$B
batch_hawkes$Pi
batch_hawkes$

station_clusts <- tibble(station = 1:m,clust = est_Z_hawkes)

stations %>% left_join(station_clusts, 
                       by = c("stat_id" = "station")) %>%
  ggplot(aes(stat_lon,
             stat_lat,color=as.factor(clust))) + 
  geom_point()

nyc_map_data <- stations %>% left_join(station_clusts, 
                                       by = c("stat_id" = "station"))


ggmap(nyc_map) + geom_point(data = nyc_map_data,
                            aes(stat_lon,                                                    stat_lat,color=as.factor(clust))) +
  theme_nothing()



clust_data <- nyc_map_data %>% filter(clust == 2)
revgeo(longitude = clust_data$stat_lon,
       latitude = clust_data$stat_lat,
       provider = 'photon',output = 'frame') 

#### should investigate the locations here which always come up ###
citi_data %>% filter(`start station longitude` %in% clust_data$stat_lon) %>%
  group_by(`start station id`) %>% tally() %>%
  arrange(desc(n))

citi_data %>%
  group_by(`start station id`) %>% tally() %>%
  arrange(desc(n))

# these don't seem to correspond to locations with more trips
# starting or ending there anyway


#### Inhomogeneous Poisson ####
window = 1/7
K <- 3 
H <- 7
dT = 1 # 
MuA_start = array(0,c(K,K,H))
tau_start = matrix(0,m,K)
system.time(non_hom_pois_est <- 
              nonhomoPois_estimator(as.matrix(bike_data),
                                    A_test,m,K,H,window,
                                    T = Time,dT,
                                    gravity = 0.001,
                                    MuA_start,tau_start) )

non_hom_pois_est$Pi
est_non_hom_pois <- apply(non_hom_pois_est$tau,1,which.max)
table(est_non_hom_pois)
station_clusts <- tibble(station = 1:m,clust = est_non_hom_pois)

stations %>% left_join(station_clusts, 
                       by = c("stat_id" = "station")) %>%
  ggplot(aes(stat_lon,
             stat_lat,color=as.factor(clust))) + 
  geom_point()


system.time(batch <- 
              batch_nonhomoPois_estimator(as.matrix(bike_data),
                                          A_test,m,K,H,window,
                                          T = Time,dT,
                                          gravity = 0.001,MuA_start,tau_start,
                                          itermax = 100,stop_eps = 0.01 ))
batch$Pi


est_Z_batch_nhp <- apply(batch$tau,1,which.max)
station_clusts <- tibble(station = 1:m,clust = est_Z_batch_nhp)

stations %>% left_join(station_clusts, 
                       by = c("stat_id" = "station")) %>%
  ggplot(aes(stat_lon,
             stat_lat,color=as.factor(clust))) + 
  geom_point()


nyc_map_data <- stations %>% left_join(station_clusts, 
                                       by = c("stat_id" = "station"))

clust_data <- nyc_map_data %>% filter(clust == 1)
revgeo(longitude = clust_data$stat_lon,
       latitude = clust_data$stat_lat,
       provider = 'photon',output = 'frame')


#### Inhom Hawkes
B_start = matrix(0,K,K)
nonHawkes_online <- nonhomoHak_estimator_eff_revised(as.matrix(bike_data),
                                                     A_test,m,K,H,window,
                                                     T = Time,dT=1,
                                                     lam = 1.75,
                                                     gravity = 0.001,
                                                     B_start,
                                                     MuA_start,tau_start,is_elbo = T)

nonHawkes_online$Pi
nonHawkes_online$elbo[length(nonHawkes_online$elbo)]*dim(bike_data)[1]


#### try restrict to Manhattan ####
zip_codes <- read_csv("C:/Users/owenw/Google Drive/Tian/Oral/data/zip_borough.csv")

manhattan <- zip_codes %>% filter(borough == "Manhattan")
manhattan

stations %>% head()

zip_only <- revgeo(longitude = stations$stat_lon,
                   latitude = stations$stat_lat,
                   provider = 'photon',output = 'frame') %>%
  select(zip)

#stations$zip <- station_zip$zip

as.numeric(as.character(station_zip$zip))

stations <- stations %>%
  mutate(zip = as.numeric(as.character(zip)))

stations %>% left_join(zip_codes) %>% select(borough)

station_zip

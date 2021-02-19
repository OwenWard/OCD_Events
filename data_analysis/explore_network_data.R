library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("cpp_files/onlineblock.cpp")

### some other data ####
fb_wall <- read_tsv(gzfile("data/facebook-wall.txt.gz"))

colnames(fb_wall) <- c("send","rec","time_stamp")

head(fb_wall)

fb_wall %>% mutate(time_stamp = lubridate::as_datetime(time_stamp)) %>%
  group_by(send,rec) %>% 
  tally() %>%
  arrange(desc(n))

### could compare this clustering with known adjacency matrix, but is that
### sort of cheating?


yt <- read_tsv(gzfile("data/youtube-links.txt.gz"))




#### Columbia Email Data ####

per_info <- read_lines(gzfile("C:/Users/owenw/Google Drive/Tian/EmailNetwork/Email-network CU/02/data0002.gz"),
                     skip_empty_rows = FALSE)

per_info <- as.data.frame(per_info)

nrow(per_info)/19

vars <- c("uni",
          "is_faculty",
          "is_student",
          "is_officer",
          "is_staff",
          "is_affiliate",
          "birth_year",
          "campus",
          "aca_field",
          "aca_department",
          "adm_department",
          "division",
          "dorm_building",
          "from_us",
          "gender",
          "postal_code",
          "school",
          "school_year",
          "student_status")


per_data <- tibble(value = per_info$per_info,variable = rep(vars,nrow(per_info)/19))

per_data %>% tail(n=20)

# this is correctly capturing everything up to here at least
unis <- per_data %>% filter(variable == "uni") %>% dplyr::select(value)

per_data$id <- rep(as.vector(unis$value),each = 19)

tidy_info <- per_data %>% pivot_wider(names_from = variable,values_from = value,
                                      values_fill = NA)

tidy_info %>% 
  select(-id) %>%
  filter(is_faculty == 1 | is_officer == 1 | is_staff == 1) %>%
  group_by(adm_department) %>% 
  count() %>% 
  arrange(desc(n))


faculty_staff <- tidy_info %>%
  select(-id) %>%
  filter(is_faculty == 1 | is_officer == 1 | is_staff == 1) %>%
  select(uni,adm_department)# %>% # could have a category for those with
  ## no dept also
  # filter(adm_department != "")


faculty_staff
##### extracting emails, using previous code ####
data_path <- "C:/Users/owenw/Google Drive/Tian/EmailNetwork/Email-network CU/02/"

files <- list.files(path = data_path)

files
# class 1 has 55 between 8 before
#days = files[5:126]  # sem 1
days = files[5:150]   # for sem 2

days_index <- c(1:146) # to account for one missing day
days_index <- days_index[c(-13)] # missing for this semester, no data


output_total <- c()


for(i in days_index){
  print(i) 
  day_sample <- readLines(gzfile(paste(data_path,days[i],sep = "")))
  output <- str_split(day_sample,pattern = ":", simplify = T)
  rec_num <- str_count(output[,3], pattern = ",")
  # only keep emails with <5 recipients
  temp <- output[which(rec_num<5),]
  # then only keep senders who are faculty...
  send <- which(temp[,2] %in% faculty_staff$uni)
  # get recipients also
  out <- str_split(temp[,3],pattern = ",")
  n_obs <- 1 # max number of recipients considered.
  seq.max <- seq_len(max(n_obs))
  tmp = t(sapply(out, "[", i = seq.max))
  if(n_obs == 1){
    tmp = t(tmp)
  }
  indices = row(tmp)
  rec <- indices[which(tmp %in% faculty_staff$uni)]
  keep <- intersect(send,rec)
  output_total <- rbind(output_total, temp[keep,])
  
}

beepr::beep()

# then only keep the first recipient of these emails

all_recips <- output_total[,3]
rec_all <- stringr::str_split(all_recips,pattern = ",",simplify = TRUE)
first_recp <- rec_all[,1]

events <- tibble(send = output_total[,2],recip = first_recp,
                 Time = output_total[,1])

events %>% 
  group_by(send,recip) %>%
  count() %>%
  arrange(desc(n))
  

summary_events <- events %>% 
  left_join(faculty_staff, by = c("send" = "uni")) %>%
  left_join(faculty_staff, by = c("recip" = "uni")) %>%
  group_by(send,recip) %>%
  count() %>%
  filter(n > 3)

summary_events

processed_events %>% ggplot(aes(Time_since)) + geom_histogram()

(max(as.numeric(events$Time))-min(as.numeric(events$Time)))/147
# corresponds to 1440 minutes in a day

start <- min(as.numeric(events$Time))

# this is in days now
events %>% mutate(Time = as.numeric(Time) - start-1,
                  Time = Time/1440 ) %>%
  ggplot(aes(Time)) + geom_histogram()

processed_events <- events %>% 
  mutate(Time = as.numeric(Time)- (start - 10) ) %>%
  mutate(Time = Time/1440) %>%
  arrange(Time)

### this is in days now, can use it almost now ###


users <- unique(c(processed_events$send,processed_events$recip))
m <- length(users)



emails <- processed_events %>% 
  mutate(Send = as.numeric(factor(send,levels = users))-1) %>% 
  mutate(Rec = as.numeric(factor(recip,levels = users))-1) %>%
  select(Send,Rec,Time)


id_code <- faculty_staff %>% mutate(new_id = as.numeric(factor(uni,levels = users))-1) %>%
  drop_na()
# this gives the same number of users as m which is good


### emails now in the format for algorithm

A = list()

for(i in 1:m){
  A[[i+1]] = i
}

for(i in 1:nrow(emails)){
  j = emails$Send[i]
  k = emails$Rec[i]
  #print(j)
  #A[[j+1]] = j
  #A[[k+1]] = k
  A[[j+1]] = c(A[[j+1]],emails$Rec[i])
}
A_test = lapply(A,unique)


# now can fit this, Poisson ####

Time <- 147
dT <- 0.5
K <- 200
Pi <- rep(1/K,K)
B <- matrix(runif(K*K),K,K)
tau <- matrix(runif(m*K),nrow=m,ncol=K)
tau <- tau/rowSums(tau)
S <- matrix(0,nrow = m,ncol = K)


results <- estimate_Poisson(full_data = as.matrix(emails),
                 A_test,m,K,Time,dT,
                 step_size = 0.5,
                 B,
                 tau,Pi,S,inter_T = 10,is_elbo = FALSE)

est_z <- apply(results$tau,1,which.max)


est_clust <- tibble(new_id = 1:m, clust = est_z) %>%
  mutate(new_id = new_id - 1) %>%
  right_join(id_code, by = "new_id")

est_clust %>% group_by(clust,adm_department) %>%
  count() %>%
  arrange(desc(n))


aricode::NMI(est_clust$clust,est_clust$adm_department)


### Inhom Hawkes ####

K <- 20
H <- 4
window <- 1/H

MuA = array(0,dim=c(K,K,H))

results_nhawkes_sim <- nonhomoHak_estimator_eff(as.matrix(emails),A_test,
                                                m,K,H,
                                                window,
                                                T=145,dT=3,
                                                lam = 0.1, gravity = 0.01,
                                                B,MuA,tau)

est_z = apply(results_nhawkes_sim$tau,1,which.max)

est_clust <- tibble(new_id = 1:m, clust = est_z) %>%
  mutate(new_id = new_id - 1) %>%
  right_join(id_code, by = "new_id")

est_clust %>% group_by(clust,adm_department) %>%
  count() %>%
  arrange(desc(n))


aricode::NMI(est_clust$clust,est_clust$adm_department)

library(here)
source(here("functions/", "utils.R"))



#### debugging the inhomogeneous Poisson model ####

## at any time point its a linear combination of the corresponding
## basis functions

t_start <- 5
t_end <- 30
window <- 5
H <- 2
period <- H * window


tvec <- rep(0, H)
h1 <- floor(t_start/window);
h2 <- floor(t_end/window);

vals <- c(h1:h2)
for(w in vals){
  cat("w ----- \n")
  print(w)
  h <- w %% H
  cat("h ----- \n")
  print(h)
  if(h1 == h2) {
    print("In here")
    tvec[h + 1] <- t_end - t_start
  }
  if(w == h1){
    print("Now in here")
    tvec[h + 1] <- tvec[h + 1] + (w + 1) * window - t_start
    print(tvec[h])
    cat("==== \n")
  }
  else if(w == h2){
    print("actually in here")
    tvec[h + 1] <- tvec[h + 1] + t_end - w * window
  }
  else{
    print("JK, in here")
    tvec[h + 1] = tvec[h + 1] + window
    print(tvec[h + 1])
    cat("==== \n")
  }
}

tvec

## whats the relationship between tvec and the realised intensity function?
## tvec a derivative?
## tvec corresponds to the compensator term,
## just sums up the time in each of the terms in the basis
m <- n
H <- 2
K <- 2
MuA <- array(runif(K * K * H), c(K, K, H))
window <- 30
tau <- matrix(1/K, nrow = m, ncol = K)
dT <- 1


results_online_inpois <- nonhomoPois_estimator(alltimes = events,
                                                         A = A_test,
                                                         m,
                                                         K,
                                                         H,
                                                         window,
                                                         Time,
                                                         dT,
                                                         gravity = 0.001,
                                                         MuA, tau)

z_true <- apply(dynppsbm$z, 2, which.max)
(ari_in_pois <- aricode::ARI(apply(results_online_inpois$tau, 1, which.max),
                             z_true))

results_online_inpois$MuA
## then try back these out to get some intensities


#### Exp 1 ####

### see if we recover the true intensity in some very simple cases

intens <- matrix(c(0.2, 1, 0.5, 0.75, 0.75, 0.5, 1, 0.25),
                 nrow = 4, byrow = TRUE)


## here just changes once during the time period

dynppsbm <- generateDynppsbmConst(intens,
                                  Time,
                                  n,
                                  prop.groups,
                                  directed = TRUE)

###
# hist(dynppsbm$data$time.seq)
proc_sim <- format_sims(sim_data = dynppsbm,
                        n = n,
                        directed = TRUE)


A_test <- proc_sim$edge
events <- proc_sim$events

H <- 2
MuA <- array(runif(K * K * H), c(K, K, H))
window <- 50 #50 ## period is H * window
dT <- 0.5
tau <- matrix(1/K, nrow = m, ncol = K)
results_online_inpois <- nonhomoPois_estimator(alltimes = events,
                                                         A = A_test,
                                                         m,
                                                         K,
                                                         H,
                                                         window,
                                                         Time,
                                                         dT,
                                                         gravity = 0.001,
                                                         MuA, tau)

z_true <- apply(dynppsbm$z, 2, which.max)
(ari_in_pois <- aricode::ARI(apply(results_online_inpois$tau, 1, which.max),
                             z_true))

results_online_inpois$MuA

results_online_inpois$MuA[1, 1, ]
results_online_inpois$MuA[1, 2, ]
results_online_inpois$MuA[2, 1, ]
results_online_inpois$MuA[2, 2, ]

intens
## so when you set the window to match the truth it does recover them 
## here, as would be expected

## how do we specify window then so that it matches these in general?


## if change the periodicity what happens
intens2 <- cbind(intens, intens)

dynppsbm2 <- generateDynppsbmConst(intens2,
                                  Time,
                                  n,
                                  prop.groups,
                                  directed = TRUE)

###
# hist(dynppsbm$data$time.seq)
proc_sim2 <- format_sims(sim_data = dynppsbm2,
                         n = n,
                         directed = TRUE)


A_test2 <- proc_sim2$edge
events2 <- proc_sim2$events

H <- 2
MuA <- array(runif(K * K * H), c(K, K, H))
window <- 25 #50 ## period is H * window
dT <- 0.5
tau <- matrix(1/K, nrow = m, ncol = K)
results_online_inpois2 <- nonhomoPois_estimator(alltimes = events2,
                                                A = A_test2,
                                                m,
                                                K,
                                                H,
                                                window,
                                                Time,
                                                dT,
                                                gravity = 0.001,
                                                MuA, tau)

z_true2 <- apply(dynppsbm2$z, 2, which.max)
(ari_in_pois2 <- aricode::ARI(apply(results_online_inpois2$tau,
                                    1, which.max),
                             z_true2))

results_online_inpois2$MuA

results_online_inpois2$MuA[1, 1, ]
results_online_inpois2$MuA[1, 2, ]
results_online_inpois2$MuA[2, 1, ]
results_online_inpois2$MuA[2, 2, ]

intens2

## so does recover them now, because window matches how often they
## change



#### Exp 3 ####

## try infer the estimated intensities
results_online_inpois$MuA
window



#' Given parameter estimates, infer the inhomogeneous Poisson intensities
#'
#' @param MuA 
#' @param window 
#' @param Time 
#'
#' @return a matrix with K*K rows and M columns, with columns 
#' corresponding to the number of jumps in the pairwise intensity function
#' @export
#'
#' @examples
infer_intensities <- function(MuA, window, t_start, t_end) {
  
  H <- dim(MuA)[3]
  period <- window * H
  
  K <- dim(MuA)[2]
  num_intens <- K * K
  h1 <- floor(t_start/window) # start closest and to the left 
  h2 <- floor(t_end/window)
  jumps <- h1:h2
  rates <- matrix(NA, nrow = K * K, ncol = length(jumps) )
  ## extra column here corresponds to new intensity at start of next
  ## window?
  
  ## then populate this piecewise rate matrix
   ## will actually give all of them at once
  row <- 1
  for(k1 in 1:K) {
    for(k2 in 1:K) {
      curr_int <- MuA[k1, k2, jumps %% H + 1]
      rates[row, ] <- curr_int
      row <- row + 1
    }
  }
  rates
}


### then want a function for plotting these also

curr_rates <- infer_intensities(MuA = results_online_inpois$MuA,
                                window = 50, t_start = 0,
                                t_end = Time)

est_rates <- curr_rates %>% 
  as_tibble() %>% 
  mutate(pair = paste0("Pair_", 1:4)) %>% 
  pivot_longer(cols = starts_with("V"), 
               ## make this general
               names_to = "Window", values_to = "rate") %>% 
  ## need to be a bit more careful with this here
  mutate(Window = as.numeric(str_extract(Window, "[:digit:]"))) %>% 
  mutate(x = t_start + (Window - 1)* window)


est_rates %>% 
  ggplot(aes(x, rate)) +
  geom_step() +
  facet_wrap(~pair, scales = "free") +
  ylim(c(0, 2)) +
  labs(y = "Intensity", x = "Time") +
  NULL

h1 <- floor(t_start/window) # start closest and to the left 
h2 <- floor(t_end/window)
jumps <- h1:h2

jumps %% H + 1

## to match the scaling of the estimated rates

true_intens <- intens[, jumps %% H + 1] %>% 
  as_tibble() %>% 
  mutate(pair = paste0("Pair_", 1:nrow(intens))) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "Window",
               values_to = "rate") %>% 
  mutate(Window = as.numeric(str_extract(Window, "[:digit:]"))) %>% 
  mutate(x = t_start + (Window - 1) * window) 

true_intens %>% 
  ggplot(aes(x, rate)) +
  geom_step() +
  facet_wrap(~pair, scales = "free") +
  xlim(c(0, Time)) +
  ylim(c(0, 1.5)) +
  labs(y = "True Intensity", x = "Time")
  


### then plot the two of them at the same time

est_rates %>% 
  ggplot(aes(x, rate)) +
  geom_step(alpha = 0.5) +
  facet_wrap(~pair) +
  geom_step(data = true_intens, aes(x, rate), colour = "Red", 
            alpha = 0.5)

## these are flipped which causes the difference



### fit to the VAST_Cell Data
library(lubridate)

data_path <- here("data/VASTchallenge08-20080315-Deinosuchus/CELL CALLS/")

calls <- read_csv(here(data_path, "CellPhoneCallRecords.csv"))

start_time <- ymd_hms("2006-06-01 00:00:00", tz = "UTC")
calls %>% 
  mutate(Time = Datetime - start_time ) %>% 
  mutate(Time_mins = as.numeric(Time)) 


### then tidy this up to the right format for our algorithm

## construct A


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
## tvec corresponds to the compensator term

H <- 3
MuA <- array(runif(K * K * H), c(K, K, H))
window <- 30
tau <- matrix(1/K, nrow = m, ncol = K)
invisible(results_online_inpois <- nonhomoPois_estimator(alltimes = events,
                                                         A = A_test,
                                                         m,
                                                         K,
                                                         H,
                                                         window,
                                                         Time,
                                                         dT,
                                                         gravity = 0.001,
                                                         MuA, tau))

(ari_in_pois <- aricode::ARI(apply(results_online_inpois$tau, 1, which.max),
                             z_true))

results_online_inpois$MuA
## then try back these out to get some intensities



### see if we recover in some simpler cases

intens <- matrix(c(0.2, 1, 0.5, 0.5, 0.5, 0.5, 1, 0.25),
                 nrow = 4, byrow = TRUE)

# intens[2, 3] <- 0.5
# intens[3, 2] <- 1

dynppsbm <- generateDynppsbmConst(intens,
                                  Time,
                                  n,
                                  prop.groups,
                                  directed = TRUE)

###
# hist(dynppsbm$data$time.seq)
proc_sim <- format_sims(sim_data = dynppsbm, n = n,
                        directed = TRUE)



### 
proc_sim <- format_sims(sim_data = dynppsbm, n = n,
                        directed = TRUE)
A_test <- proc_sim$edge
events <- proc_sim$events

H <- 2
MuA <- array(runif(K * K * H), c(K, K, H))
window <- 50
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

## try infer the estimated intensities


results_online_inpois$MuA
window



#' Given parameter estimates, infer the inhomogeneous Poisson intensities
#'
#' @param MuA 
#' @param window 
#' @param Time 
#'
#' @return
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
}


### then want a function for plotting these also

tmp <- rates %>% 
  as_tibble() %>% 
  mutate(pair = paste0("Pair_", 1:4)) %>% 
  pivot_longer(cols = V1:V3, names_to = "Window", values_to = "rate") %>% 
  ## need to be a bit more careful with this here
  mutate(Window = as.numeric(str_extract(tmp$Window, "[:digit:]")))


tmp %>% 
  mutate(x = t_start + (Window - 1)* window) %>% 
  ggplot(aes(x, rate)) +
  geom_step() +
  facet_wrap(~pair, scales = "free")

tmp %>% 
  ggplot(aes(Window, rate)) +
  geom_step() +
  facet_wrap(~pair)

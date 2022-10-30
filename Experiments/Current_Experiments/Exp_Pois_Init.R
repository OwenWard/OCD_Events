##### 27th October, 2022
### Aim of these experiments is to assess the impact 
### of the initialization procedures as we make our data
### "sparse", in terms of the rates of the underlying point processes

library(here)

source(here("functions/", "utils.R"))
source(here("functions/init_fcn.R"))

### given data from an inhomogeneous Poisson, estimate the 
### corresponding parameters

## here will need the window size and the number of basis functions H
## look at Guanhua's cpp code for this also which should be useful


window <- 10
H <- 2

period <- window * H


## simulate some data, fit the model to it

avec <- c(1, 5)
t_start <- 0
t_end <- 100
window <- 10

sim_data <- genpois_nonhomo(t_start, T = t_end, window, avec)

as.vector(sim_data)

hist(sim_data)

# then want to fit to this, given window, H, accounting for 
# repeated windows, etc

num_in_wind <- function(events, jumps, which_h, window){
  ## given the breakpoints (here right is upper limit,)
  ## return vector of length length(jumps)-1,
  ## giving how many events in a specific window 
  counts <- rep(0, length(jumps) - 1)
  time_H <- rep(0,H)
  counts_H <- rep(0, H)
  for(i in seq_along(jumps[-length(jumps)])){
    lower <- jumps[i]
    upper <- jumps[i+1]
    counts[i] <- length(which(events < upper & events > lower))
    time_H[which_h[i]] <- time_H[which_h[i]] + window 
  }
  list(counts = counts, time = time_H)
}

fit_inpois <- function(t_start, t_end, events, window, H){
  
  ## bin the data based on windows, then just concatenate
  ## it and fit separately?
  ## this won't work for hawkes because history is common,
  ## or would that just be an additional step then?
  
  num_events <- rep(0, H)
  # obs_time <- rep(0, H)
  h1 <- floor(t_start/window) # start closest and to the left 
  h2 <- floor(t_end/window) # this will throw out incomplete intervals
  jumps <- h1:h2 * window # to get the actual jumps
  which_h <- rep(1:H, length.out = length(jumps) - 1)
  ## then bin the data based on these jumps
  count_data <- num_in_wind(events, jumps, which_h, window)
  obs_time <- 
  
}

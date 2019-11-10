# times down to the minute
time_data <- readRDS("C:/Users/owenw/Documents/Online_Point_Process/Sem6tidy.RDS")
all_users <- unique(c(time_data$Sender_id,time_data$Rec_id))

m = length(all_users)

# fit the model?
# construct A
##
A = list()



for(i in 1:nrow(time_data)){
  j = time_data$Sender_id[i]
  k = time_data$Rec_id[i]
  #print(j)
  A[[j+1]] = j
  A[[k+1]] = k
  A[[j+1]] = c(A[[j+1]],time_data$Rec_id[i])
}


A_test = lapply(A,unique)

dim(time_data)
#B[[3568]]




library(Rcpp)
library(RcppArmadillo)


sourceCpp("C:/Users/owenw/Documents/Online_Point_Process/onlineblock.cpp")



### Poisson Simulation ###
T = 90
dT = 0.5
K <- 3
Mu <- matrix(c(0.6,0.2,0.3,0.1,1.0,0.4,0.5,0.4,0.8),K,K,byrow = TRUE)
B <- matrix(0,K,K,byrow = TRUE)
m <- 100
Pi <- matrix(c(0.4,0.3,0.3),1,3)
Z <- c(rep(0,m*Pi[1]),rep(1,m*Pi[2]),rep(2,m*Pi[3]))

system.time(alltimes <- sampleBlockHak(T, A, Z, Mu, B, lam = 1))


###
T = 90
dT = 0.5
m = length(all_users)
Pi = c(0.3,0.3,0.4)
B = matrix(c(1.2,0.5,0.5,0.5,1.1,.65,0.75,0.85,1.15),nrow = K,ncol = K,byrow = T)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)



results <- estimate_Poisson(full_data = as.matrix(time_data),tau,B,Pi,S,A_test,m,K,dT,T)
results$B

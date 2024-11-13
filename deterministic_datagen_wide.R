library(doParallel)
library(foreach)

# Calculate the number of cores
getDoParWorkers()                                          
detectCores()                                                      
cl=makeCluster(7)                                              
registerDoParallel(cl)                                                               
getDoParWorkers()   

myfunc = function(m)
{
  library(MASS);library(ResourceSelection);
  library(dplyr); library(glm2);
  library(data.table)
  setDTthreads(1)
  
  source("datagen.R")
  set.seed(1129)
  seeds = floor(runif(1000)*10^8);
  set.seed(seeds[m])
  
rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
N <- 100000
theta0=-1; theta1=-1; theta2=1; theta3=1; theta4 = -2;
sigma=1

W = runif(N, 0, 1)
L =  rbinom(N, 1, 0.5)
A = rbinom(N, 3, plogis(-2+L+W))
B = rbinom(N, 3, plogis(1-2*L+W))
C = rbinom(N, 3, plogis(-1+L+W))

Aint = rep(0,N)
Bint = ifelse(B + A > 3, 3, B + A)
Y = rbinom(N,1, plogis(theta0+theta1*Aint+theta2*Bint+theta3*C+theta4*L))

fit = glm2(Y~L, family=binomial())
myparam  = c(coef(fit))

return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"true.csv")

stopCluster(cl)
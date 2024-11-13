library(doParallel)
library(foreach)

# Calculate the number of cores
getDoParWorkers()                                          
detectCores()                                                      
cl=makeCluster(8)                                              
registerDoParallel(cl)                                                               
getDoParWorkers()   

#for(m in 1:sim)
myfunc = function(m)
{
  library(MASS);library(ResourceSelection);
  library(dplyr); library(glm2);
  library(data.table)
  setDTthreads(1)
  
  logit <- function(term) {
    return( ifelse(!is.na(term),log(term/(1-term)),NA) )
  }
  
  EXPIT <- function(term) {
    return( ifelse(!is.na(term),exp(term)/(1+exp(term)),NA) )
  }
  
  source("datagen.R")
  set.seed(1129)
  seeds = floor(runif(1000)*10^8);
  set.seed(seeds[m])
  
  n <- 5000
  theta0=-1; theta1=-1; theta2=1; theta3=1; theta4 = -2;
  sigma=1
  df <- datagen(N=n, sigma=sigma,
            theta0=theta0, theta1=theta1, theta2=theta2, 
            theta3=theta3, theta4=theta4)
  
  tmpdata = df

  ##################
  ######time 2######
  ##################
  yfit = glm2(Y ~ A+B+C+L, family=binomial(), data = tmpdata) ; 
  ydat = tmpdata; ydat$B=ifelse(ydat$A+ydat$B>3, 3, ydat$A+ydat$B); ydat$A=0; 
  ydat$ypred = predict(yfit, newdata = ydat, type="response")

  fitmsm = glm2(ypred ~ L, family=quasibinomial(), data = ydat) ; 
  myparam = coef(fitmsm)
  myparam
  
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"seq_reg.csv")

stopCluster(cl)
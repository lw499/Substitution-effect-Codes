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
  library(data.table); library(nnet);
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
  
  bfit =  multinom(B ~ as.factor(A) * as.factor(C) * L, data = df) 
  afit = multinom(A ~ as.factor(C) * L, data=df)
  
  df$pred_diff = predict(bfit, newdata = df, type="probs")
  df$denom1 = ifelse(df$B==0, df$pred_diff[,1], df$pred_diff[,4])
  df$denom1 = ifelse(df$B==1, df$pred_diff[,2], df$denom1)
  df$denom1 = ifelse(df$B==2, df$pred_diff[,3], df$denom1)

  df$pred_diffa = predict(afit, newdata = df, type="probs")
  df$denom2 = ifelse(df$A==0, df$pred_diffa[,1], df$pred_diffa[,4])
  df$denom2 = ifelse(df$A==1, df$pred_diffa[,2], df$denom2)
  df$denom2 = ifelse(df$A==2, df$pred_diffa[,3], df$denom2)

  ##a=0
  df0 = df; df0$A=0;
  df$pred_diff = predict(bfit, newdata = df0, type="probs")
  diff0 = df$B-0
  num1_0 = ifelse(diff0==0, df$pred_diff[,1], 0)
  num1_0 = ifelse(diff0==1, df$pred_diff[,2], num1_0)
  num1_0 = ifelse(diff0==2, df$pred_diff[,3], num1_0)
  num1_0 = ifelse(diff0==3, df$pred_diff[,4], num1_0)

  cumsum = 1
  num2_0 = I(df$B==3)*(1-cumsum)
  num_0 = (num1_0 + num2_0)*df$pred_diffa[,1]
    
  ##a=1
  df1 = df; df1$A=1;
  df$pred_diff = predict(bfit, newdata = df1, type="probs")
  diff1 = df$B-1
  num1_1 = ifelse(diff1==0, df$pred_diff[,1], 0)
  num1_1 = ifelse(diff1==1, df$pred_diff[,2], num1_1)
  num1_1 = ifelse(diff1==2, df$pred_diff[,3], num1_1)

  cumsum = df$pred_diff[,1] + df$pred_diff[,2] + df$pred_diff[,3] 
  num2_1 = I(df$B==3)*(1-cumsum)
  num_1 = (num1_1 + num2_1)*df$pred_diffa[,2]
  
  ##a=2
  df2 = df; df2$A=2;
  df$pred_diff = predict(bfit, newdata = df2, type="probs")
  diff2 = df$B-2
  num1_2 = ifelse(diff2==0, df$pred_diff[,1], 0)
  num1_2 = ifelse(diff2==1, df$pred_diff[,2], num1_2)

  cumsum = df$pred_diff[,1] + df$pred_diff[,2]
  num2_2 = I(df$B==3)*(1-cumsum)
  num_2 = (num1_2 + num2_2)*df$pred_diffa[,3]
  
  ##a=3
  df3 = df; df3$A=3;
  df$pred_diff = predict(bfit, newdata = df3, type="probs")
  diff3 = df$B-3
  num1_3 = ifelse(diff3==0, df$pred_diff[,1], 0)

  cumsum = df$pred_diff[,1]
  num2_3 = I(df$B==3)*(1-cumsum)
  num_3 = (num1_3 + num2_3)*df$pred_diffa[,4]
  
  ## total
  df$num_tot = num_0 + num_1 + num_2 + num_3
  
  ##################
  ######result######
  ##################
  
  test = glm2(Y~L, weight = (I(A==0)*num_tot)/(denom1*denom2), family=quasibinomial(), data=df)
  myparam = coef(test)
  myparam
  
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"ipw.csv")

stopCluster(cl)
library(doParallel)
library(foreach)

# Calculate the number of cores
getDoParWorkers()                                          
detectCores()                                                      
cl=makeCluster(30)                                              
registerDoParallel(cl)                                                               
getDoParWorkers()   

#for(m in 1:sim)
myfunc = function(m)
{
  library(MASS);library(ResourceSelection);
  library(dplyr); library(glm2); library(hal9001);
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
  df <- datagen(N=n, sigma=sigma)
  
  df$A0 = ifelse(df$A==0,1,0); df$B0 = ifelse(df$B==0,1,0)
  df$A1 = ifelse(df$A==1,1,0); df$B1 = ifelse(df$B==1,1,0)
  df$A2 = ifelse(df$A==2,1,0); df$B2 = ifelse(df$B==2,1,0)
  df$A3 = ifelse(df$A==3,1,0); df$B3 = ifelse(df$B==3,1,0)
  
  xmatB0 = as.matrix(cbind(df$L1, df$L2,  df$A, df$C)); 
  fitB0 = fit_hal(X = xmatB0, Y = df$B0, family = "binomial")
  dftmp = df[df$B>0,]; xmatB1 = as.matrix(cbind(dftmp$L1, dftmp$L2, dftmp$A, dftmp$C));
  fitB1 = fit_hal(X = xmatB1, Y = dftmp$B1, family = "binomial")
  dftmp = df[df$B>1,]; xmatB2 = as.matrix(cbind(dftmp$L1, dftmp$L2, dftmp$A, dftmp$C));
  fitB2 = fit_hal(X = xmatB2, Y = dftmp$B2, family = "binomial")
    
  xmatB = as.matrix(cbind(df$L1, df$L2,  df$A, df$C)); 
  df$Bpred0 = predict(fitB0, new_data=xmatB); predicted0 = df$Bpred0;
  df$Bpred1 = predict(fitB1, new_data=xmatB); predicted1 = df$Bpred1*(1-df$Bpred0);
  df$Bpred2 = predict(fitB2, new_data=xmatB); predicted2 = df$Bpred2*(1-df$Bpred1)*(1-df$Bpred0);
  predicted3 = 1-predicted0-predicted1-predicted2
  
  xmatA0 = as.matrix(cbind(df$L1, df$L2,  df$C)); 
  fitA0 = fit_hal(X = xmatA0, Y = df$A0, family = "binomial")
  dftmp = df[df$B>0,]; xmatA1 = as.matrix(cbind(dftmp$L1, dftmp$L2, dftmp$A, dftmp$C));
  fitA1 = fit_hal(X = xmatA1, Y = dftmp$A1, family = "binomial")
  dftmp = df[df$B>1,]; xmatA2 = as.matrix(cbind(dftmp$L1, dftmp$L2, dftmp$A, dftmp$C));
  fitA2 = fit_hal(X = xmatA2, Y = dftmp$A2, family = "binomial")
  
  xmatA = as.matrix(cbind(df$L1, df$L2,  df$C))
  df$Apred0 = predict(fitA0, new_data=xmatA); predicted0a = df$Apred0;
  df$Apred1 = predict(fitA1, new_data=xmatA); predicted1a = df$Apred1*(1-df$Apred0);
  df$Apred2 = predict(fitA2, new_data=xmatA); predicted2a = df$Apred2*(1-df$Apred1)*(1-df$Apred0);
  predicted3a = 1-predicted0a-predicted1a-predicted2a
  
  df$denom1 = ifelse(df$B==0, predicted0, predicted3)
  df$denom1 = ifelse(df$B==1, predicted1, df$denom1)
  df$denom1 = ifelse(df$B==2, predicted2, df$denom1)
  
  df$denom2 = ifelse(df$A==0, predicted0a, predicted3a)
  df$denom2 = ifelse(df$A==1, predicted1a, df$denom2)
  df$denom2 = ifelse(df$A==2, predicted2a, df$denom2)
  
  ##a=0
  df0 = df; df0$A=0; xmatB = as.matrix(cbind(df0$L1, df0$L2, df0$A, df0$C))
  df$Bpred0 = predict(fitB0, new_data=xmatB); predictdiff0 = df$Bpred0;
  df$Bpred1 = predict(fitB1, new_data=xmatB); predictdiff1 = df$Bpred1*(1-df$Bpred0);
  df$Bpred2 = predict(fitB2, new_data=xmatB); predictdiff2 = df$Bpred2*(1-df$Bpred1)*(1-df$Bpred0);
  predictdiff3 = 1-predictdiff0-predictdiff1-predictdiff2
  
  diff0 = df$B-0
  num1_0 = ifelse(diff0==0, predictdiff0, 0)
  num1_0 = ifelse(diff0==1, predictdiff1, num1_0)
  num1_0 = ifelse(diff0==2, predictdiff2, num1_0)
  num1_0 = ifelse(diff0==3, predictdiff3, num1_0)
  
  cumsum = 1
  num2_0 = I(df$B==3)*(1-cumsum)
  num_0 = (num1_0 + num2_0)*predicted0a
  
  ##a=1
  df1 = df; df1$A=1; xmatB = as.matrix(cbind(df1$L1,df1$L2, df1$A, df1$C))
  df$Bpred0 = predict(fitB0, new_data=xmatB); predictdiff0 = df$Bpred0;
  df$Bpred1 = predict(fitB1, new_data=xmatB); predictdiff1 = df$Bpred1*(1-df$Bpred0);
  df$Bpred2 = predict(fitB2, new_data=xmatB); predictdiff2 = df$Bpred2*(1-df$Bpred1)*(1-df$Bpred0);
  predictdiff3 = 1-predictdiff0-predictdiff1-predictdiff2
  
  diff1 = df$B-1
  num1_1 = ifelse(diff1==0, predictdiff0, 0)
  num1_1 = ifelse(diff1==1, predictdiff1, num1_1)
  num1_1 = ifelse(diff1==2, predictdiff2, num1_1)
  
  cumsum = predictdiff0 + predictdiff1 + predictdiff2
  num2_1 = I(df$B==3)*(1-cumsum)
  num_1 = (num1_1 + num2_1)*predicted1a
  
  ##a=2
  df2 = df; df2$A=2; xmatB = as.matrix(cbind(df2$L1, df2$L2, df2$A, df2$C))
  df$Bpred0 = predict(fitB0, new_data=xmatB); predictdiff0 = df$Bpred0;
  df$Bpred1 = predict(fitB1, new_data=xmatB); predictdiff1 = df$Bpred1*(1-df$Bpred0);
  df$Bpred2 = predict(fitB2, new_data=xmatB); predictdiff2 = df$Bpred2*(1-df$Bpred1)*(1-df$Bpred0);
  predictdiff3 = 1-predictdiff0-predictdiff1-predictdiff2
  
  diff2 = df$B-2
  num1_2 = ifelse(diff2==0, predictdiff0, 0)
  num1_2 = ifelse(diff2==1, predictdiff1, num1_2)
  
  cumsum = predictdiff0 + predictdiff1
  num2_2 = I(df$B==3)*(1-cumsum)
  num_2 = (num1_2 + num2_2)*predicted2a
  
  ##a=3
  df3 = df; df3$A=3; xmatB = as.matrix(cbind(df3$L1, df3$L2, df3$A, df3$C))
  df$Bpred0 = predict(fitB0, new_data=xmatB); predictdiff0 = df$Bpred0;
  df$Bpred1 = predict(fitB1, new_data=xmatB); predictdiff1 = df$Bpred1*(1-df$Bpred0);
  df$Bpred2 = predict(fitB2, new_data=xmatB); predictdiff2 = df$Bpred2*(1-df$Bpred1)*(1-df$Bpred0);
  predictdiff3 = 1-predictdiff0-predictdiff1-predictdiff2
  
  diff3 = df$B-3
  num1_3 = ifelse(diff3==0, predictdiff0, 0)
  
  cumsum = predictdiff0
  num2_3 = I(df$B==3)*(1-cumsum)
  num_3 = (num1_3 + num2_3)*predicted3a
  
  ## total
  df$num_tot = num_0 + num_1 + num_2 + num_3
  tmpdata = df
  
  ##################
  ######time 2######
  ##################
  ydat = tmpdata
  xmatY = as.matrix(cbind(ydat$L1, ydat$L2, ydat$A, ydat$B, ydat$C)); 
  yfito = fit_hal(X = xmatY, Y = ydat$Y, family = "binomial")
  ydat$ypredo = qlogis(predict(yfito, new_data = xmatY))
  yfitup = glm(Y~L2, offset = ypredo, weight = (I(A==0)*num_tot)/(denom1*denom2), family=quasibinomial(), data=ydat)
  
  ydat = tmpdata; ydat$B=ifelse(ydat$A+ydat$B>3, 3, ydat$A+ydat$B); ydat$A=0; 
  xmatYup = as.matrix(cbind(ydat$L1, ydat$L2, ydat$A, ydat$B, ydat$C)); 
  ydat$ypredo = qlogis(predict(yfito, new_data = xmatYup))
  ydat$ypred = predict(yfitup, newdata = ydat, type="response")

  fitmsm = glm2(ypred ~ L2, family=quasibinomial(), data = ydat) ; 
  myparam = coef(fitmsm)
  myparam
  
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"dr_ML.csv")

stopCluster(cl)
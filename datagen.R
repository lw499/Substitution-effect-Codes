#time varying case:
#L1: private insurance
#L2: screening of rectal STI (yes/no)
#L3: log transformed CD4 (maximal value of 10)

datagen <- function(N, sigma,
                    theta0=theta0, theta1=theta1, theta2=theta2, 
                    theta3=theta3, theta4=theta4)
{
  id <- seq(1, N, by=1) # Define individual id number
  # Generate baseline common cause of time-varying covariates
  W = runif(N, 0, 1)
  L =  rbinom(N, 1, 0.5)
  A = rbinom(N, 3, plogis(-2+L+W))
  B = rbinom(N, 3, plogis(1-2*L+W))
  C = rbinom(N, 3, plogis(-1+L+W))
  
  Y = rbinom(N,1, plogis(theta0+theta1*A+theta2*B+theta3*C+theta4*L))
  
  # Consolidate data in a single data frame
  temp_data <- data.frame(id, L, A, B, C, Y)
  return(temp_data)
}



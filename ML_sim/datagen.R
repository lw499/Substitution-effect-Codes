#time varying case:
#L1: private insurance
#L2: screening of rectal STI (yes/no)
#L3: log transformed CD4 (maximal value of 10)

datagen <- function(N, sigma)
{
  theta0=-1; theta1=-1; theta2=1; theta3=1; theta4 = -2; theta5 = 1; theta6=-0.25; theta7 = 0.25;
  id <- seq(1, N, by=1) # Define individual id number
  # Generate baseline common cause of time-varying covariates
  W = runif(N, 0, 1)
  L1 =  rbinom(N, 1, 0.5)
  L2 = rnorm(N, 0, 1)
  A = rbinom(N, 3, plogis(-2+0.75*L1 + 0.5*L2+W))
  B = rbinom(N, 3, plogis(-1+0.5*L1 - 0.2*exp(0.5*L2)+W))
  C = rbinom(N, 3, plogis(-1+L1+W))
  Y = rbinom(N,1, plogis(theta0+theta1*A+theta2*B+theta3*C+theta4*L1+theta5*L2+theta6*L2^2+theta7*L1*L2))
  
  # Consolidate data in a single data frame
  temp_data <- data.frame(id, L1, L2, A, B, C, Y)
  return(temp_data)
}



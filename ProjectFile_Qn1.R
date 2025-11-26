# Statistical Inference

library(ggplot2)

# ==============================================================================
# 1 - Exponential data: estimation bias and efficiency
# ==============================================================================


#-------------------------------------------------------------------------------
# i)
#-------------------------------------------------------------------------------
# Using the function rexp, create a function that generates a sample of exponential
# data and computes the maximum likelihood estimator (MLE) of λ. Call this function
# mle exp. The inputs of this function should be the value of λ and the sample size, n.
# Roll No: 25189506,25203533

# MLE  λ = n/ΣXi
# CRLB = (λ^2)/n

mle_exp <- function(lambda, sample_size) {
  exp_data <- rexp(sample_size, rate = lambda)
  lambda_hat <- 1/mean(exp_data)
  return(lambda_hat)
}


#-------------------------------------------------------------------------------
# ii)
#-------------------------------------------------------------------------------
# Generate one MLE using mle exp with λ = 1 and n = 10. Note that any time you run
# this, you will get different results. Therefore, to ensure that I (the examiner)
# can replicate the exact result which appears in your report, insert 
# set.seed(SEEDNUMBER) in your code just prior to mle exp where SEEDNUMBER is the 
# sum of your two student ID numbers.
# For λ = 1
seednumber <- 25189506+25203533
set.seed(seednumber)
lambda <- 1
sample_size <- 10
lambda_hat <- mle_exp(lambda, sample_size)
lambda_hat



#-------------------------------------------------------------------------------
# iii)
#-------------------------------------------------------------------------------
# The result of (ii) only applies to one sample so does not tell us much about
# the general properties of the estimator; for this we need replicate samples.
# Using a for loop, repeat part (ii) 1000 times. Now calculate the mean and
# variance of these 1000 MLEs. Also calculate the mean squared error. Note: just
# prior to your loop, insert set.seed(SEEDNUMBER)
seednumber <- 25189506+25203533
set.seed(seednumber)
sum_lambda_hat = 0
sum_lambda_hat_2 = 0
#lambda_hat_list <- list()
for (x in 1:1000) {
  lambda_hat <- mle_exp(lambda, sample_size)
  sum_lambda_hat <-  sum_lambda_hat + lambda_hat
  sum_lambda_hat_2 <-  sum_lambda_hat_2 + (lambda_hat**2)
  #lambda_hat_list <- append(lambda_hat_list, lambda_hat)
}
E_lambda <- sum_lambda_hat/1000
E_lambda2 <- sum_lambda_hat_2/1000
Variance_lambda <- E_lambda2 - (E_lambda**2)
#Var <- var(unlist(lambda_hat_list))
Bias <- E_lambda-lambda
MSE <- Variance_lambda + (Bias)**2


#-------------------------------------------------------------------------------
# iv)
#-------------------------------------------------------------------------------
# Repeat step (iii) but for n = 20, n = 50, n = 200, n = 400, n = 1000 (again
# insert set.seed(SEEDNUMBER) just prior to each loop),and produce the following
# graphs:
# a) bias versus sample size,
# b) efficiency versus sample size (note: you need to calculate the CRLB
#    by hand first), and
# c) mean squared error versus sample size.
# In each of these graphs, indicate the values of the analytic bias, efficiency, 
# and mean squared error, i.e., those computed by hand. Show the hand calculations
# for these quantities and also for the CRLB.

lambda <- 1
SampleSizeList <- list(20,50,200,400,1000)
BiasList <- list()
MseList <- list()
EfficiencyList <- list()
for (sample_size in SampleSizeList) {
  seednumber <- 25189506+25203533
  set.seed(seednumber)
  sum_lambda_hat = 0
  sum_lambda_hat_2 = 0
  for (x in 1:1000) {
    lambda_hat <- mle_exp(lambda, sample_size)
    sum_lambda_hat <-  sum_lambda_hat + lambda_hat
    sum_lambda_hat_2 <-  sum_lambda_hat_2 + (lambda_hat**2)
  }
  E_lambda <- sum_lambda_hat/1000
  E_lambda2 <- sum_lambda_hat_2/1000
  Variance_lambda <- E_lambda2 - (E_lambda**2)
  Bias <- E_lambda-lambda
  MSE <- Variance_lambda + (Bias)**2
  BiasList <- append(BiasList,Bias)
  MseList <- append(MseList,MSE)
  CRLB = (lambda**2) / sample_size
  Efficiency = CRLB/Variance_lambda
  EfficiencyList <- append(EfficiencyList, Efficiency)
}

print(BiasList)
print(MseList)

# Analytical Values - Bias, Variance, MSE, Efficiency
Analytical_Bias_List <- list()
Analytical_Variance_List <- list()
Analytical_MSE_List <- list()
Analytical_Efficiency_List <- list()
for (sample_size in SampleSizeList) {
  bias <- lambda/(sample_size-1)
  Analytical_Bias_List <- append(Analytical_Bias_List,bias)
  variance <- (lambda**2 * sample_size**2)/((sample_size-1)**2 * (sample_size-2))
  Analytical_Variance_List <- append(Analytical_Variance_List,variance)
  MSE <- variance + (bias**2)
  Analytical_MSE_List <- append(Analytical_MSE_List,MSE)
  CRLB = (lambda**2) / sample_size
  Efficiency <- CRLB/variance
  Analytical_Efficiency_List <- append(Analytical_Efficiency_List,Efficiency)
}


# Final Data
OutputData <- list(
  SampleSize = unlist(SampleSizeList),
  Bias = unlist(BiasList),
  AnalyticalBias = unlist(Analytical_Bias_List),
  MSE = unlist(MseList),
  AnalyticalMSE = unlist(Analytical_MSE_List),
  Efficiency = unlist(EfficiencyList),
  AnalyticalEfficiency = unlist(Analytical_Efficiency_List)
)

df <- as.data.frame(OutputData)

ggplot(data=df, aes(x=SampleSize, y=Bias)) + geom_line() + geom_point()
ggplot(data=df, aes(x=SampleSize, y=MSE)) + geom_line() + geom_point()
ggplot(data=df, aes(x=SampleSize, y=Efficiency)) + geom_line() + geom_point()

ggplot() +
  geom_line(data = df, aes(SampleSize, Bias), color = "red") +
  geom_line(data = df, aes(SampleSize, AnalyticalBias), color = "blue")
  
ggplot() +
  geom_line(data = df, aes(SampleSize, MSE), color = "red") +
  geom_line(data = df, aes(SampleSize, AnalyticalMSE), color = "blue")

ggplot() +
  geom_line(data = df, aes(SampleSize, Efficiency), color = "red") +
  geom_line(data = df, aes(SampleSize, AnalyticalEfficiency), color = "blue")

#-------------------------------------------------------------------------------
# v)
#-------------------------------------------------------------------------------
# Now repeat step (iv) but for λ = 0.5 and λ = 2, and do not forget 
# set.seed(SEEDNUMBER). Discuss the findings for the various λ values and sample
# sizes investigated.

# For λ = 0.5
lambda <- 0.5
SampleSizeList <- list(20,50,200,400,1000)
BiasList <- list()
MseList <- list()
EfficiencyList <- list()
for (sample_size in SampleSizeList) {
  seednumber <- 25189506+25203533
  set.seed(seednumber)
  sum_lambda_hat = 0
  sum_lambda_hat_2 = 0
  for (x in 1:1000) {
    lambda_hat <- mle_exp(lambda, sample_size)
    sum_lambda_hat <-  sum_lambda_hat + lambda_hat
    sum_lambda_hat_2 <-  sum_lambda_hat_2 + (lambda_hat**2)
  }
  E_lambda <- sum_lambda_hat/1000
  E_lambda2 <- sum_lambda_hat_2/1000
  Variance_lambda <- E_lambda2 - (E_lambda**2)
  Bias <- E_lambda-lambda
  MSE <- Variance_lambda + (Bias)**2
  BiasList <- append(BiasList,Bias)
  MseList <- append(MseList,MSE)
  CRLB = (lambda**2) / sample_size
  Efficiency = CRLB/Variance_lambda
  EfficiencyList <- append(EfficiencyList, Efficiency)
}

print(BiasList)
print(MseList)

OutputData <- list(
  SampleSize = unlist(SampleSizeList),
  Bias = unlist(BiasList),
  MSE = unlist(MseList),
  Efficiency = unlist((EfficiencyList))
)

df <- as.data.frame(OutputData)

ggplot(data=df, aes(x=SampleSize, y=Bias)) + geom_point()
ggplot(data=df, aes(x=SampleSize, y=MSE)) + geom_point()
ggplot(data=df, aes(x=SampleSize, y=Efficiency)) + geom_point()


# For λ = 2 
lambda <- 2
SampleSizeList <- list(20,50,200,400,1000)
BiasList <- list()
MseList <- list()
EfficiencyList <- list()
for (sample_size in SampleSizeList) {
  seednumber <- 25189506+25203533
  set.seed(seednumber)
  sum_lambda_hat = 0
  sum_lambda_hat_2 = 0
  for (x in 1:1000) {
    lambda_hat <- mle_exp(lambda, sample_size)
    sum_lambda_hat <-  sum_lambda_hat + lambda_hat
    sum_lambda_hat_2 <-  sum_lambda_hat_2 + (lambda_hat**2)
  }
  E_lambda <- sum_lambda_hat/1000
  E_lambda2 <- sum_lambda_hat_2/1000
  Variance_lambda <- E_lambda2 - (E_lambda**2)
  Bias <- E_lambda-lambda
  MSE <- Variance_lambda + (Bias)**2
  BiasList <- append(BiasList,Bias)
  MseList <- append(MseList,MSE)
  CRLB = (lambda**2) / sample_size
  Efficiency = CRLB/Variance_lambda
  EfficiencyList <- append(EfficiencyList, Efficiency)
}

print(BiasList)
print(MseList)

OutputData <- list(
  SampleSize = unlist(SampleSizeList),
  Bias = unlist(BiasList),
  MSE = unlist(MseList),
  Efficiency = unlist((EfficiencyList))
)

df <- as.data.frame(OutputData)

ggplot(data=df, aes(x=SampleSize, y=Bias)) + geom_point() + geom_smooth()
ggplot(data=df, aes(x=SampleSize, y=MSE)) + geom_point()  + geom_smooth()
ggplot(data=df, aes(x=SampleSize, y=Efficiency)) + geom_point()  + geom_smooth()


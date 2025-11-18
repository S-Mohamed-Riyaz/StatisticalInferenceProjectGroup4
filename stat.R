  #1.1
mle_exp <- function(lambda, n) {
  x <- rexp(n, rate = lambda)
  return(1 / mean(x))
}



#1.2
SEEDNUMBER <-25203533 + 25189506
set.seed(SEEDNUMBER)
one_mle <- mle_exp(lambda = 1, n = 10)
one_mle


#1.3
set.seed(SEEDNUMBER)
mles <- numeric(1000)
for (i in 1:1000) {
  mles[i] <- mle_exp(lambda=1, n=10)
}



mean_mles <- mean(mles)
var_mles <- var(mles)
mse_mles  <- mean( (mles - 1)^2 )  

mean_mles
var_mles
mse_mles
#1.4
sample_sizes <- c(20, 50, 200, 400, 1000)

results <- data.frame(
  n   = sample_sizes,
  mean = NA,
  var  = NA,
  mse  = NA,
  bias = NA,
  efficiency = NA
)


CRLB <- function(lambda, n) lambda^2 / n

for (i in 1:length(sample_sizes)) {
  n <- sample_sizes[i]
  
  set.seed(SEEDNUMBER)
  mles <- numeric(1000)
  for (j in 1:1000) {
    x <- rexp(n, rate = 1)
    mles[j] <- 1 / mean(x)   
  }
  
  
  results$mean[i] <- mean(mles)
  results$var[i]  <- var(mles)
  results$bias[i] <- mean(mles) - 1
  results$mse[i]  <- mean( (mles - 1)^2 )
  results$efficiency[i] <- CRLB(1, n) / var(mles)
}

results

library(ggplot2)
library(tidyr)


results_long <- results %>%
  pivot_longer(cols = c("bias", "efficiency", "mse"),
               names_to = "metric",
               values_to = "value")


ggplot(results_long, aes(x = n, y = value)) +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    title = "Bias, Efficiency and MSE vs Sample Size",
    x = "Sample size n",
    y = ""
  ) +
  theme_bw(base_size = 14)

#1.5


compute_results <- function(lambda_value, sample_sizes, SEEDNUMBER = 25203533) {
    
  CRLB <- function(lambda, n) lambda^2 / n
  
  results <- data.frame(
    n   = sample_sizes,
    mean = NA,
    var  = NA,
    mse  = NA,
    bias = NA,
    efficiency = NA
  )
  
  for (i in 1:length(sample_sizes)) {
    n <- sample_sizes[i]
    
    set.seed(SEEDNUMBER)
    mles <- numeric(1000)
    for (j in 1:1000) {
      x <- rexp(n, rate = lambda_value)
      mles[j] <- 1 / mean(x)
    }
    
    results$mean[i] <- mean(mles)
    results$var[i]  <- var(mles)
    results$bias[i] <- mean(mles) - lambda_value
    results$mse[i]  <- mean((mles - lambda_value)^2)
    results$efficiency[i] <- CRLB(lambda_value, n) / var(mles)
  }
  
  return(results)
}



#lambda = 0.5
results_lambda05 <- compute_results(lambda_value = 0.5,
                                    sample_sizes = sample_sizes,
                                    SEEDNUMBER = SEEDNUMBER)

results_lambda05



#lambda = 2 
results_lambda2 <- compute_results(lambda_value = 2,
                                   sample_sizes = sample_sizes,
                                   SEEDNUMBER = SEEDNUMBER)

results_lambda2


results_long_05 <- results_lambda05 %>%
  pivot_longer(cols = c("bias", "efficiency", "mse"),
               names_to = "metric",
               values_to = "value")

ggplot(results_long_05, aes(x = n, y = value)) +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    title = "Bias, Efficiency and MSE vs Sample Size (lambda = 0.5)",
    x = "Sample size n",
    y = ""
  ) +
  theme_bw(base_size = 14)


results_long_2 <- results_lambda2 %>%
  pivot_longer(cols = c("bias", "efficiency", "mse"),
               names_to = "metric",
               values_to = "value")

ggplot(results_long_2, aes(x = n, y = value)) +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    title = "Bias, Efficiency and MSE vs Sample Size (lambda = 2)",
    x = "Sample size n",
    y = ""
  ) +
  theme_bw(base_size = 14)




# ==============================================================================
# 2- Comparing maximum likelihood and method of moments
# ==============================================================================

library(dplyr)
library(ggplot2)
#-------------------------------------------------------------------------------
# i)
#-------------------------------------------------------------------------------
# Consider the distribution with density function given by f(x) = cx^(θ−1)
# where θ > 0 and x ∈ [0,xmax]. Randomly select a value for xmax to be used
# throughout this question by running sample(2:15, size = 1) (with set.seed
# prior to it for reproducibility). Note that this is quite similar to Q3 in
# Tutorial 2 but xmax ̸ = 1 here.

seednumber <- 25189506+25203533
set.seed(seednumber)

xmax <- sample(2:15, size = 1)  #xmax = 8

# f(x) = θx^(θ−1) / (8^θ)

#-------------------------------------------------------------------------------
# ii)
#-------------------------------------------------------------------------------
# By hand, compute the MLE and the MoM estimator for θ. Show all workings.

# MOM θ_hat = Xbar / (8-Xbar)

# MLE θ_hat = n / (nlog8 - Σlogxi)


#-------------------------------------------------------------------------------
# iii)
#-------------------------------------------------------------------------------
# Unlike Q1, this is a non-standard distribution, which you must simulate from
# yourself. This can be done using the method of inversion described here.
# First, derive the cumulative distribution function (cdf), F(x) = x0
# f(u)du, and then the associated inverse cdf, F−1(x).
# Next, the random variable X can be generated as X = F−1(U) where
# U ∼Uniform(0,1), and uniform variables can be generated in R via the
# runif function. Using this approach, generate samples of size n = 5000
# for a few different values of θ values, and, using a histogram, describe
# the shape of the distribution and how it changes with respect to θ; also
# indicate the theoretical expected value E(X) on the histogram.

# x = 8u^(1/θ)
# EX = 8θ/(θ+1)

# θ = 0.5
theta_real <- 0.5

seednumber <- 25189506+25203533
set.seed(seednumber)

u <- runif(5000)
x <- xmax*u^(1/theta_real)
EX <- xmax*theta_real/(theta_real+1)

ggplot(data.frame(x),aes(x=x)) + geom_histogram() + 
  annotate("text", x = Inf, y = Inf,
           label = paste0("EX = ", round(EX,2)),
           hjust = 1.1, vjust = 1.5,
           size = 5)

# θ = 1
theta_real <- 1

seednumber <- 25189506+25203533
set.seed(seednumber)

u <- runif(5000)
x <- xmax*u^(1/theta_real)
EX <- xmax*theta_real/(theta_real+1)

ggplot(data.frame(x),aes(x=x)) + geom_histogram() + 
  annotate("text", x = Inf, y = Inf,
           label = paste0("EX = ", round(EX,2)),
           hjust = 1.1, vjust = 1.5,
           size = 5)

# θ = 2
theta_real <- 2

seednumber <- 25189506+25203533
set.seed(seednumber)

u <- runif(5000)
x <- xmax*u^(1/theta_real)
EX <- xmax*theta_real/(theta_real+1)

ggplot(data.frame(x),aes(x=x)) + geom_histogram() + 
  annotate("text", x = Inf, y = Inf,
           label = paste0("EX = ", round(EX,2)),
           hjust = 1.1, vjust = 1.5,
           size = 5)

# θ = 4
theta_real <- 4

seednumber <- 25189506+25203533
set.seed(seednumber)

u <- runif(5000)
x <- xmax*u^(1/theta_real)
EX <- xmax*theta_real/(theta_real+1)

ggplot(data.frame(x),aes(x=x)) + geom_histogram() + 
  annotate("text", x = Inf, y = Inf,
           label = paste0("EX = ", round(EX,2)),
           hjust = 1.1, vjust = 1.5,
           size = 5)


#-------------------------------------------------------------------------------
# iv)
#-------------------------------------------------------------------------------
# Generate one sample of size n = 100 with θ = 1 and compute both the MLE
# and MoM.

# MOM, θ_hat = Xbar / (8 - Xbar)
# MLE, θ_hat = n / (nlog8 - ΣlogXi) 

theta_real <- 1
n <- 100

seednumber <- 25189506+25203533
set.seed(seednumber)

u <- runif(n)
x <- xmax*u^(1/theta_real)

sum_logxi <- sum(log(x))

Xbar <- mean(x) 

MOM_Theta_hat <- Xbar / (xmax-Xbar)

MLE_Theta_hat <- n / ((n*log(xmax)) - sum_logxi)

#-------------------------------------------------------------------------------
# v)
#-------------------------------------------------------------------------------
# Now generate two additional samples: one of size n = 10 and one of 
# size n = 1000. For all three samples, plot the log-likelihood function with
# respect to θ, i.e., θ values on the x-axis and the ℓ(θ) values on the y-axis.
# Indicate where both the MLE and MoM lie on these plots. Describe how the MLE 
# and MoM compare to each other, and to the true value of θ, in all three cases.

# l(θ) = nlogθ + (θ-1)ΣlogXi - nθlog8

theta_real <- 1
n_values <- c(10, 100, 1000)

seednumber <- 25189506+25203533
set.seed(seednumber)

sample_data <- list()
estimator_results <- data.frame()
data_list <- list()

log_likelihood <- function(theta, X) {
  n <- length(X)
  sum_ln_X <- sum(log(X))
  
  n * log(theta) - n * theta * log(xmax) + (theta - 1) * sum_ln_X
}


# The previously generated sample (n=100) is now included in the loop
for (n in n_values) {
  u <- runif(n)
  x <- xmax*u^(1/theta_real)
  
  # Calculate Estimators
  sum_logxi <- sum(log(x))
  Xbar <- mean(x)
  
  MOM_Theta_hat <- Xbar / (xmax-Xbar)
  
  MLE_Theta_hat <- n / ((n*log(xmax)) - sum_logxi)
  
  # Finding the log likelihood values for MLE and MOM estimator
  ll_mle <- log_likelihood(MLE_Theta_hat, x)
  ll_mom <- log_likelihood(MOM_Theta_hat, x)
  
  # Store estimators for the plot
  new_row <- data.frame(
    n = n,
    MLE = MLE_Theta_hat,
    MoM = MOM_Theta_hat,
    True = theta_real,
    LL_MLE = ll_mle,
    LL_MoM = ll_mom
  )
  estimator_results <- rbind(estimator_results, new_row)
  
  # Generate log-likelihood values over a range of theta
  theta_range <- seq(0.5, xmax, length.out = 5000)
  ll_values <- log_likelihood(theta_range, x)
  
  # Store plot data
  sample_data[[as.character(n)]] <- data.frame(
    Theta = theta_range,
    LogLikelihood = ll_values,
    n_size = as.factor(n)
  )
}

plot_data <- bind_rows(sample_data)


# Plotting the data
p <- ggplot(plot_data, aes(x = Theta, y = LogLikelihood)) +
  geom_line(aes(group = n_size), color = "darkblue", size = 0.8) +
  
  # Add MLE points (Y-value from LL_MLE column)
  geom_point(data = estimator_results, 
             aes(x = MLE, y = LL_MLE, color = "MLE"), 
             size = 3, shape = 17) + # Triangle shape
  
  # Add MoM points (Y-value from LL_MoM column)
  geom_point(data = estimator_results, 
             aes(x = MoM, y = LL_MoM, color = "MoM"), 
             size = 3, shape = 18) + # Diamond shape
  
  # Add True Theta line
  geom_vline(aes(xintercept = True, linetype = "True θ"), 
             data = distinct(estimator_results, n, True), 
             color = "black", size = 0.5) +
  
  # Separate plots by sample size
  facet_wrap(~ n, ncol = 3, scales = "free_y", 
             labeller = labeller(n = 
                                   c("10" = "Sample Size n = 10",
                                     "100" = "Sample Size n = 100",
                                     "1000" = "Sample Size n = 1000"))) +
  
  # Customize labels and theme
  labs(
    title = expression(paste("Log-Likelihood Function ", l(theta), " for Different Sample Sizes (True ", theta, " = 1)")),
    x = expression(theta),
    y = expression(l(theta)),
    color = "Estimator",
    linetype = ""
  ) +
  scale_color_manual(values = c("MLE" = "red", "MoM" = "orange")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p)


#-------------------------------------------------------------------------------
# vi)
#-------------------------------------------------------------------------------
# For each of the three log-likelihood plots, draw a horizontal line that is
# exactly d/2 units below the maximum where d = Χ2(1,0.05). Drop two vertical
# lines down from where the horizontal hits the log-likelihood function to form
# a Wilks confidence interval for θ

# d =  Χ2(1,0.05) = 3.841
# d/2 = 1.9205

d = 3.841
interval = d/2 

theta_real <- 1
n_values <- c(10, 100, 1000)

seednumber <- 25189506+25203533
set.seed(seednumber)

sample_data <- list()
estimator_results <- data.frame()
data_list <- list()

log_likelihood <- function(theta, X) {
  n <- length(X)
  sum_ln_X <- sum(log(X))
  
  n * log(theta) - n * theta * log(xmax) + (theta - 1) * sum_ln_X
}


# The previously generated sample (n=100) is now included in the loop
for (n in n_values) {
  u <- runif(n)
  x <- xmax*u^(1/theta_real)
  
  # Calculate Estimator
  sum_logxi <- sum(log(x))
  
  MLE_Theta_hat <- n / ((n*log(xmax)) - sum_logxi)
  
  # Finding the log likelihood values for MLE estimator
  ll_mle <- log_likelihood(MLE_Theta_hat, x)
  
  # Calculate Wilks CI threshold
  ll_threshold <- ll_mle - interval
  
  # Store estimators for the plot
  new_row <- data.frame(
    n = n,
    MLE = MLE_Theta_hat,
    True = theta_real,
    LL_MLE = ll_mle,
    LL_Threshold = ll_threshold
  )
  estimator_results <- rbind(estimator_results, new_row)
  
  # Generate log-likelihood values over a range of theta
  theta_range <- seq(0.5, xmax, length.out = 5000)
  ll_values <- log_likelihood(theta_range, x)
  
  # Store plot data
  sample_data[[as.character(n)]] <- data.frame(
    Theta = theta_range,
    LogLikelihood = ll_values,
    n_size = as.factor(n)
  )
}

plot_data <- bind_rows(sample_data)


# Plotting the data
p <- ggplot(plot_data, aes(x = Theta, y = LogLikelihood)) +
  geom_line(aes(group = n_size), color = "darkblue", size = 0.8) +
  
  # Add MLE points (Y-value from LL_MLE column)
  geom_point(data = estimator_results, 
             aes(x = MLE, y = LL_MLE, color = "MLE"), 
             size = 3, shape = 17) + # Triangle shape + 
  
  # Add True Theta line
  geom_vline(aes(xintercept = True, linetype = "True θ"), 
             data = distinct(estimator_results, n, True), 
             color = "black", size = 0.5) +
  
  geom_hline(aes(yintercept = LL_Threshold, color = "95% CI Threshold"), 
             data = distinct(estimator_results, n, LL_Threshold), 
             linetype = "dotted", size = 1) +
  
  # Separate plots by sample size
  facet_wrap(~ n, ncol = 3, scales = "free_y", 
             labeller = labeller(n = 
                                   c("10" = "Sample Size n = 10",
                                     "100" = "Sample Size n = 100",
                                     "1000" = "Sample Size n = 1000"))) +
  
  # Customize labels and theme
  labs(
    title = expression(paste("Log-Likelihood Function ", l(theta), " for Different Sample Sizes (True ", theta, " = 1)")),
    x = expression(theta),
    y = expression(l(theta)),
    color = "Estimator",
    linetype = ""
  ) +
  scale_color_manual(values = c("MLE" = "red")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p)


#-------------------------------------------------------------------------------
# vii)
#-------------------------------------------------------------------------------
# For each of the three sample sizes, find the MLE numerically using the numerical
# optimiser nlm. Note: this is done by minimising minus the log-likelihood (which
# is equivalent to maximising the log-likelihood). When using nlm, set 
# hessian = TRUE so that nlm stores the so-called “hessian”, i.e., the second 
# derivative of the objective function. Since we minimise minus the log-likelihood,
# this is already the observed information. Use this to produce a 95% Wald 
# confidence interval for θ. How do these compare with the Wilks confidence 
# intervals from part (v)?

theta_real <- 1
n_values <- c(10, 100, 1000)

seednumber <- 25189506+25203533
set.seed(seednumber)

ci_results <- data.frame()
data_list <- list()

# Z-score for 95% CI (two-tailed) (For Walds CI)
z_95percent <- qnorm(0.975)

# d =  Χ2(1,0.05) = 3.841 (For Wilks CI)
d = 3.841
interval = d/2

log_likelihood <- function(theta, X) {
  n <- length(X)
  sum_ln_X <- sum(log(X))
  
  n * log(theta) - n * theta * log(xmax) + (theta - 1) * sum_ln_X
}

negative_log_likelihood <- function(theta, X) {
  -log_likelihood(theta,X)
}

for (n in n_values) {
  u <- runif(n)
  x <- xmax*u^(1/theta_real)
  
  # Calculate Estimator
  sum_logxi <- sum(log(x))
  MLE_Theta_hat <- n / ((n*log(xmax)) - sum_logxi)
  
  theta_assumed = MLE_Theta_hat
  
  # Run nlm to minimize -l(θ) and to get the observed Information I(θ).
  nlm_result <- nlm(
    f = negative_log_likelihood, 
    p = theta_assumed, 
    X = x,            
    hessian = TRUE
  )
  
  # Finding the numerical MLE
  Numerical_theta_hat <- nlm_result$estimate
  
  # Getting the observed information I(θ)
  I_theta <- nlm_result$hessian[1, 1]
  
  mean <- theta_assumed
  variance <- 1/I_theta
  sd <- variance**0.5
  
  # Walds Confidence Interval
  walds_upperlimit <- mean + z_95percent*sd
  walds_lowerlimit <- mean - z_95percent*sd
  
  # Finding the log likelihood values for MLE estimator
  ll_mle <- log_likelihood(MLE_Theta_hat, x)
  
  # Calculate Wilks CI threshold
  ll_threshold <- ll_mle - interval
  
  # Defining a function for the equation to be solved
  equation_fn <- function(theta, X, ll_threshold) {
    return(log_likelihood(theta, X) - ll_threshold)
  }
  
  # Finding the roots of the equations
  wilks_lowerlimit <- uniroot(
    f = equation_fn, 
    interval = c(0.01, MLE_Theta_hat), 
    X = x, 
    ll_threshold = ll_threshold
  )$root
  
  wilks_upperlimit <- uniroot(
    f = equation_fn, 
    interval = c(MLE_Theta_hat,8), 
    X = x, 
    ll_threshold = ll_threshold
  )$root
  
  # Store estimators for the plot
  new_row <- data.frame(
    n = n,
    MLE_Analytical = MLE_Theta_hat,
    MLE_Numerical = Numerical_theta_hat,
    True = theta,
    Walds_CI_upper = walds_upperlimit,
    Walds_CI_lower = walds_lowerlimit,
    Wilks_CI_upper = wilks_upperlimit,
    Wilks_CI_lower = wilks_lowerlimit
  )
  ci_results <- rbind(ci_results, new_row)
}

print(ci_results)


#-------------------------------------------------------------------------------
# viii)
#-------------------------------------------------------------------------------
# Repeat step (vii) 1000 times. In other words, use nlm to numerically maximise
# the log-likelihood function and compute 95% Wald confidence intervals in each
# of 1000 samples at each of the three sample sizes (but there is no need to
# compute Wilks intervals). Also compute the MoM. How do the MLEs compare to the
# true θ value at each of the three sample sizes? And how do the MoM estimators
# compare? Check the coverage of the Wald confidence intervals for each of the
# three sample sizes, i.e., what proportion of the intervals contain the
# true θ value? Are the results as expected?

theta_real <- 1
n_values <- c(10, 100, 1000)

seednumber <- 25189506+25203533
set.seed(seednumber)

# Z-score for 95% CI (two-tailed) (For Walds CI)
z_95percent <- qnorm(0.975)

overall_results <- data.frame()
coverage_results <- data.frame()

for (n in n_values) {
  sample_results <- data.frame(
    MLE_Theta = numeric(1000),
    MoM_Theta = numeric(1000),
    Wald_LowerLimit = numeric(1000),
    Wald_UpperLimit = numeric(1000)
  )
  
  for (i in 1:1000){
    u <- runif(n)
    x <- xmax*u^(1/theta_real)
    
    # Calculate MOM Estimator
    Xbar <- mean(x)
    MOM_Theta_hat <- Xbar / (xmax-Xbar)
    sample_results$MOM_Theta[i] <- MOM_Theta_hat
    
    # Calculate MLE Estimator
    sum_logxi <- sum(log(x))
    MLE_Theta_hat <- n / ((n*log(xmax)) - sum_logxi)
    sample_results$MLE_Theta[i] <- MLE_Theta_hat
    
    theta_assumed = MLE_Theta_hat
    
    # Run nlm to minimize -l(θ) and to get the observed Information I(θ).
    nlm_result <- nlm(
      f = negative_log_likelihood, 
      p = theta_assumed, 
      X = x,            
      hessian = TRUE
    )
    
    # Finding the numerical MLE
    Numerical_theta_hat <- nlm_result$estimate
    
    # Getting the observed information I(θ)
    I_theta <- nlm_result$hessian[1, 1]
    
    mean <- theta_assumed
    variance <- 1/I_theta
    sd <- variance**0.5
    
    # Walds Confidence Interval
    walds_upperlimit <- mean + z_95percent*sd
    walds_lowerlimit <- mean - z_95percent*sd
    
    sample_results$Wald_UpperLimit[i] <- walds_upperlimit
    sample_results$Wald_LowerLimit[i] <- walds_lowerlimit
  }
  sample_results$n <- as.factor(n)
  overall_results <- bind_rows(overall_results, sample_results)
  
}


results_summary <- overall_results %>% group_by(n) %>%
  summarise(
    Mean_MLE = mean(MLE_Theta),
    Mean_MOM = mean(MOM_Theta),
    True = theta_true,
    Coverage = mean(Wald_LowerLimit <= theta_true & Wald_UpperLimit >= theta_true)
  )

print(results_summary)





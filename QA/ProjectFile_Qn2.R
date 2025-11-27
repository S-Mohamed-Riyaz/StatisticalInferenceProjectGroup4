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

ggplot(data.frame(x), aes(x)) +
  geom_histogram(fill = "lightblue", color = "white", bins = 30) +
  geom_vline(aes(xintercept = EX),
             color = "red",
             linewidth = 1,
             linetype = "dashed") +
  annotate("label",
           x = max(x), y = Inf,
           label = paste0("E[X] = ", round(EX, 2)),
           hjust = 1.1, vjust = 1.2,
           size = 4, fill = "white") +
  theme_minimal(base_size = 14)

# θ = 1
theta_real <- 1

seednumber <- 25189506+25203533
set.seed(seednumber)

u <- runif(5000)
x <- xmax*u^(1/theta_real)
EX <- xmax*theta_real/(theta_real+1)

ggplot(data.frame(x), aes(x)) +
  geom_histogram(fill = "lightblue", color = "white", bins = 30) +
  geom_vline(aes(xintercept = EX),
             color = "red",
             linewidth = 1,
             linetype = "dashed") +
  annotate("label",
           x = max(x), y = Inf,
           label = paste0("E[X] = ", round(EX, 2)),
           hjust = 1.1, vjust = 1.2,
           size = 4, fill = "white") +
  theme_minimal(base_size = 14)

# θ = 2
theta_real <- 2

seednumber <- 25189506+25203533
set.seed(seednumber)

u <- runif(5000)
x <- xmax*u^(1/theta_real)
EX <- xmax*theta_real/(theta_real+1)

ggplot(data.frame(x), aes(x)) +
  geom_histogram(fill = "lightblue", color = "white", bins = 30) +
  geom_vline(aes(xintercept = EX),
             color = "red",
             linewidth = 1,
             linetype = "dashed") +
  annotate("label",
           x = max(x), y = Inf,
           label = paste0("E[X] = ", round(EX, 2)),
           hjust = 1.1, vjust = 1.2,
           size = 4, fill = "white") +
  theme_minimal(base_size = 14)

# θ = 4
theta_real <- 4

seednumber <- 25189506+25203533
set.seed(seednumber)

u <- runif(5000)
x <- xmax*u^(1/theta_real)
EX <- xmax*theta_real/(theta_real+1)

ggplot(data.frame(x), aes(x)) +
  geom_histogram(fill = "lightblue", color = "white", bins = 30) +
  geom_vline(aes(xintercept = EX),
             color = "red",
             linewidth = 1,
             linetype = "dashed") +
  annotate("label",
           x = max(x), y = Inf,
           label = paste0("E[X] = ", round(EX, 2)),
           hjust = 1.1, vjust = 1.2,
           size = 4, fill = "white") +
  theme_minimal(base_size = 14)


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
MOM_Theta_hat
MLE_Theta_hat <- n / ((n*log(xmax)) - sum_logxi)
MLE_Theta_hat
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
  theta_range <- seq(0.0001, xmax, length.out = 5000)
  ll_values <- log_likelihood(theta_range, x)
  
  # Store plot data
  sample_data[[as.character(n)]] <- data.frame(
    Theta = theta_range,
    LogLikelihood = ll_values,
    n_size = as.factor(n)
  )
}

# Plotting the data
plot_data <- bind_rows(sample_data)

plot_data$n_size <- factor(plot_data$n_size, levels = as.character(n_values))

# Join estimator results so each panel has only its own points
plot_data <- plot_data %>%
  left_join(estimator_results %>% 
              rename(n_size = n) %>% 
              mutate(n_size = as.character(n_size)),
            by = "n_size")

# Plot
p <- ggplot(plot_data, aes(x = Theta, y = LogLikelihood)) +
  
  # Likelihood curve for this sample size
  geom_line(color = "darkblue", size = 0.8) +
  
  # MLE point for this sample size only
  geom_point(aes(x = MLE, y = LL_MLE, color = "MLE"),
             size = 3, shape = 17, na.rm = TRUE) +
  
  # MoM point for this sample size only
  geom_point(aes(x = MoM, y = LL_MoM, color = "MoM"),
             size = 3, shape = 18, na.rm = TRUE) +
  
  # True θ line for this sample size only
  geom_vline(aes(xintercept = True),
             linetype = "dashed", color = "black", size = 0.5, na.rm = TRUE) +
  
  # Separate panels for each sample size
  facet_wrap(~ n_size, ncol = 3, scales = "free_y",
             labeller = labeller(n_size = 
                                   c("10" = "Sample Size n = 10",
                                     "100" = "Sample Size n = 100",
                                     "1000" = "Sample Size n = 1000"))) +
  
  labs(
    title = expression(paste("Log-Likelihood Function ", l(theta), " for Different Sample Sizes (True ", theta, " = 1)")),
    x = expression(theta),
    y = expression(l(theta)),
    color = "Estimator"
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

d = qchisq(0.95, df = 1)

# Prepare a list to store CI bounds for each n
ci_list <- list()

for (n in n_values) {
  df <- sample_data[[as.character(n)]]
  est <- estimator_results %>% filter(n == !!n)
  
  # Maximum log-likelihood at MLE
  LL_max <- est$LL_MLE
  theta_hat <- est$MLE
  
  # Horizontal line
  hline <- LL_max - d/2
  
  # Split likelihood curve
  df_left  <- df %>% filter(Theta <= theta_hat)
  df_right <- df %>% filter(Theta >= theta_hat)
  
  # Interpolate the points where LL = LL_max - d/2
  ci_lower <- approx(df_left$LogLikelihood, df_left$Theta, xout = hline)$y
  ci_upper <- approx(df_right$LogLikelihood, df_right$Theta, xout = hline)$y
  
  ci_list[[as.character(n)]] <- data.frame(
    n_size = as.character(n),
    hline = hline,
    lower = ci_lower,
    upper = ci_upper
  )
}

ci_df <- bind_rows(ci_list)
ci_df$n_size <- factor(ci_df$n_size, levels = as.character(n_values))

# Combine plot data with estimator results
plot_data <- bind_rows(sample_data, .id = "n_size")
plot_data$n_size <- factor(plot_data$n_size, levels = as.character(n_values))
plot_data <- plot_data %>%
  left_join(estimator_results %>% 
              rename(n_size = n) %>% 
              mutate(n_size = as.character(n_size)),
            by = "n_size")

# Plot with Wilks CI
p <- ggplot(plot_data, aes(x = Theta, y = LogLikelihood)) +
  geom_line(color = "darkblue", size = 0.8) +
  geom_point(aes(x = MLE, y = LL_MLE, color = "MLE"), size = 3, shape = 17, na.rm = TRUE) +
  geom_vline(aes(xintercept = True), linetype = "dashed", color = "black", size = 0.5, na.rm = TRUE) +
  
  # Wilks horizontal lines
  geom_hline(data = ci_df, aes(yintercept = hline), linetype = "dotted", color = "darkgreen") +
  
  # Wilks vertical lines at CI bounds
  geom_vline(data = ci_df, aes(xintercept = lower), linetype = "dashed", color = "darkgreen") +
  geom_vline(data = ci_df, aes(xintercept = upper), linetype = "dashed", color = "darkgreen") +
  
  facet_wrap(~ n_size, ncol = 3, scales = "free_y",
             labeller = labeller(n_size = 
                                   c("10" = "Sample Size n = 10",
                                     "100" = "Sample Size n = 100",
                                     "1000" = "Sample Size n = 1000"))) +
  
  labs(
    title = expression(paste("Log-Likelihood Function and Wilks CI for Different Sample Sizes")),
    x = expression(theta),
    y = expression(l(theta)),
    color = "Estimator"
  ) +
  scale_color_manual(values = c("MLE" = "red")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p)

#-------------------------------------------------------------------------------

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
    p = theta_real, 
    X = x,            
    hessian = TRUE
  )
  
  # Finding the numerical MLE
  Numerical_theta_hat <- nlm_result$estimate
  
  # Getting the observed information I(θ)
  I_theta <- nlm_result$hessian[1, 1]
  
  mean <- Numerical_theta_hat
  variance <- 1/I_theta
  sd <- variance**0.5
  
  # Walds Confidence Interval
  walds_upperlimit <- mean + z_95percent*sd
  walds_lowerlimit <- mean - z_95percent*sd

  # Wilks Confidence Interval from part (vi)
  wilks_lowerlimit <- ci_df$lower[ci_df$n_size == n]
  wilks_upperlimit <- ci_df$upper[ci_df$n_size == n]
  
  # Store estimators for the plot
  new_row <- data.frame(
    n = n,
    MLE = Numerical_theta_hat,
    True = theta_real,
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
    True = theta_real,
    Coverage = mean(Wald_LowerLimit <= theta_real & Wald_UpperLimit >= theta_real)
  )

print(results_summary)





#3.2
SEEDNUMBER <-25203533 + 25189506
set.seed(SEEDNUMBER)  
n <- 100
lambda <- 1


x <- rexp(n, rate = lambda)



boot_medians <- numeric(1000)

for (i in 1:1000) {
  boot_sample <- sample(x, size = n, replace = TRUE)
  boot_medians[i] <- median(boot_sample)
}

boot_ci <- quantile(boot_medians, probs = c(0.025, 0.975))
boot_ci


#3.3


R <- 1000      
m_true <- log(2) 

cover <- logical(R)
lengths <- numeric(R)

for (r in 1:R) {

  x <- rexp(n, rate = lambda)
  boot_meds <- numeric(1000)
  for (i in 1:1000) {
    boot_sample <- sample(x, size = n, replace = TRUE)
    boot_meds[i] <- median(boot_sample)
  }
  
  ci <- quantile(boot_meds, probs = c(0.025, 0.975))
  L <- ci[1]; U <- ci[2]
  cover[r] <- (L <= m_true && m_true <= U)
  lengths[r] <- U - L
}

coverage_est <- mean(cover)
avg_length <- mean(lengths)
coverage_est
avg_length

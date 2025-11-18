SEEDNUMBER <-25203533
set.seed(SEEDNUMBER)
xmax <- sample(2:15, 1)

simulate_custom <- function(n, theta, xmax) {
  u <- runif(n)
  x <- xmax * u^(1/theta)
  return(x)
}


thetas <- c(0.5, 1, 2, 5)

samples <- lapply(thetas, function(th) simulate_custom(5000, th, xmax))
names(samples) <- paste0("theta_", thetas)

samples
par(mfrow = c(2,2))  # 4 graphiques

for (theta in thetas) {
  x <- simulate_custom(5000, theta, xmax)
  
  hist(x, breaks = 40, col = "lightblue",
       main = paste("theta =", theta),
       xlab = "x", probability = TRUE)
  
  # moyenne théorique
  EX <- theta / (theta + 1) * xmax
  
  abline(v = EX, col = "red", lwd = 2)
  legend("topright", legend = paste0("E[X] = ", round(EX, 2)),
         col = "red", lwd = 2)
}


library(nnet)
#---------------------------------------------------------------------iii
data <- read.csv("garments_worker_productivity.csv")
data <- na.omit(data)


Y <- data$actual_productivity
data$quarter    <- factor(data$quarter)
data$department <- factor(data$department)
data$day        <- factor(data$day)
data$team       <- factor(data$team)


hist(Y, breaks = 20, main = "Actual productivity", xlab = "actual_productivity")

one_level_cols <- names(which(sapply(data, function(col) length(unique(col))) == 1))
data <- data[, !(names(data) %in% one_level_cols)]


X <- model.matrix(actual_productivity ~ ., data = data)[, -1]
n <- nrow(X)
p <- ncol(X)
summary(data)

num_vars <- c("targeted_productivity","smv","wip","over_time",
              "incentive","idle_time","idle_men",
              "no_of_style_change","no_of_workers")

cor(data[, num_vars], use="complete.obs")
plot(data$targeted_productivity, Y, 
     xlab="targeted_productivity", ylab="actual_productivity",
     main="Relation between targeted product. and actual product.")


boxplot(Y, horizontal = TRUE, main="Boxplot of Y")
boxplot(Y ~ department, data = data, main="Productivity by department")
boxplot(Y ~ day, data = data, main="Productivity by day")
boxplot(Y ~ quarter, data = data, main="Productivity by quarter")


#------------------------------------------------------------------------------------iv
best_mse <- rep(NA, 6)
best_bic <- rep(NA, 6)

for (q in 1:6) {
  mse_list <- c()
  
  for (r in 1:5) {      # seulement 5 runs pour simplifier
    set.seed(r + q)
    fit <- nnet(X, Y, size=q, linout=TRUE, maxit=1000, trace=FALSE)
    
    pred <- predict(fit, X)
    mse_list[r] <- mean((Y - pred)^2)
  }
  
  mse_q <- min(mse_list)
  best_mse[q] <- mse_q
  
  k <- (p + 2)*q + 1
  loglik <- -n/2 * (log(mse_q) + log(2*pi) + 1)
  best_bic[q] <- -2*loglik + log(n)*(k+1)
}

results <- data.frame(q=1:6, MSE=best_mse, BIC=best_bic)
results

lm_fit <- lm(Y ~ X)
lm_pred <- predict(lm_fit)
lm_mse <- mean((Y - lm_pred)^2)

k_lm <- length(coef(lm_fit))
loglik_lm <- -n/2 * (log(lm_mse) + log(2*pi) + 1)
lm_bic <- -2*loglik_lm + log(n)*(k_lm + 1)

lm_mse
lm_bic


#------------------


q_best <- 1

set.seed(25189506+25203533)
full_model <- nnet(X, Y, size=q_best, linout=TRUE, maxit=1000, trace=FALSE)
full_pred  <- predict(full_model, X)
full_mse   <- mean((Y - full_pred)^2)

k_full <- (p + 2)*q_best + 1
loglik_full <- -n/2 * (log(full_mse) + log(2*pi) + 1)
bic_full <- -2*loglik_full + log(n)*(k_full + 1)

bic_full
#----------------------------v
var_names <- colnames(X)
bic_drop  <- rep(NA, p)

for (j in 1:p) {
  X_drop <- X[, -j, drop=FALSE]
  p_drop <- ncol(X_drop)
  
  set.seed(25189506+25203533)
  fit_j <- nnet(X_drop, Y, size=q_best, linout=TRUE, maxit=1000, trace=FALSE)
  pred_j <- predict(fit_j, X_drop)
  mse_j <- mean((Y - pred_j)^2)
  
  k_j <- (p_drop + 2)*q_best + 1
  loglik_j <- -n/2 * (log(mse_j) + log(2*pi) + 1)
  bic_drop[j] <- -2*loglik_j + log(n)*(k_j + 1)
}


importance <- data.frame(
  variable = var_names,
  BIC_without = bic_drop,
  Delta_BIC = bic_drop - bic_full
)

importance <- importance[order(-importance$Delta_BIC), ]
importance

#---------------------------------------------------------------vi
plot_effect <- function(var_name, model, X) {
  j <- which(colnames(X) == var_name)
  grid <- seq(min(X[, j]), max(X[, j]), length.out = 50)
  tau <- rep(NA, 50)
  
  for (k in 1:50) {
    Xtmp <- X
    Xtmp[, j] <- grid[k]
    tau[k] <- mean(predict(model, Xtmp))
  }
  
  plot(grid, tau, type="l", lwd=2,
       xlab=var_name, ylab="tau_hat",
       main=paste("Effet moyen de", var_name))
}
plot_effect_dummy <- function(var_name, model, X) {
  j <- which(colnames(X) == var_name)
  
  tau0 <- mean(predict(model, replace(X, j, 0)))
  tau1 <- mean(predict(model, replace(X, j, 1)))
  
  plot(c(0,1), c(tau0, tau1), type="b", lwd=2, pch=16,
       xlab=var_name, ylab="tau_hat",
       main=paste("Effet de", var_name))
}

top3 <- head(importance$variable, 3)
top3

for (v in top3) plot_effect(v, full_model, X)


#----------------------------------------------------------vii
mu_nn <- predict(full_model, X)
phi_hat <- mean((Y - mu_nn)^2)  
lower_nn <- mu_nn - 1.96 * sqrt(phi_hat)
upper_nn <- mu_nn + 1.96 * sqrt(phi_hat)

head(data.frame(Y, mu_nn, lower_nn, upper_nn))

pred_lm <- predict(lm_fit, interval="prediction", level=0.95)
head(pred_lm)


                
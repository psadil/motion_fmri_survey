library(dplyr)

ledoit_wolf_shrinkage <- function(X){
  n_samples <- nrow(X)
  n_features <- ncol(X)
  
  X <- scale(X, scale = FALSE)
  X2 <- X^2
  emp_cov_trace <- colSums(X2) / n_samples
  mu <- sum(emp_cov_trace) / n_features
  beta_ <- sum(crossprod(X2))
  delta_ <- sum(crossprod(X)^2) / n_samples^2
  # use delta_ to compute beta
  beta <- 1.0 / (n_features * n_samples) * (beta_ / n_samples - delta_)
  # delta is the sum of the squared coefficients of (<X.T,X> - mu*Id) / p
  delta <- (delta_ - 2.0 * mu * sum(emp_cov_trace) + n_features * mu^2) / n_features
  # get final beta as the min between beta and delta
  # We do this to prevent shrinking more than "1", which would invert
  # the value of covariances
  beta <- min(beta, delta)
  # finally get shrinkage
  if (beta == 0){
    shrinkage <- 0
  }else{
    shrinkage <- beta / delta
  }
  shrinkage
}


lw <- function(X){
  shrinkage <- ledoit_wolf_shrinkage(X)
  n_samples <- nrow(X)
  n_features <- ncol(X)
  Xc <- scale(X, scale = FALSE)
  # emp_cov <- cov(X)
  emp_cov <- crossprod(Xc) / n_samples
  mu <- sum(psych::tr(emp_cov)) / n_features
  shrunk_cov <- (1.0 - shrinkage) * emp_cov
  shrunk_cov[seq(1, length(emp_cov), by=n_features+1)] <- shrunk_cov[seq(1, length(emp_cov), by=n_features+1)] + shrinkage * mu
  shrunk_cov
}

n_samples <- 100
n_features <- 100
# Sigma <- (1-m)*diag(n_features) + m*U*U'


# X <- MASS::mvrnorm(100, mu=rep(0,100), Sigma = )

library(mgcv)

n.iter <- 1000

sim_fun <- function(out_misspecify = FALSE, samp_misspecify = FALSE) {
  
  X <- runif(1000, 0, 10)
  U <- exp((X - 5)/2)
  
  mu <- 1 + X + 0.5*(X^2)
  Y <- rnorm(1000, mu, sqrt(X))

  eta <- 1.5 - 0.3*X - 0.15*X^2
  S <- rbinom(1000, 1, plogis(eta))
  D <- rep(c(0,1), each = 500)
  
  if (samp_misspecify == TRUE) {
    sampmod <- glm(S ~ X, family = binomial(), subset = D == 1)
  } else {
    sampmod <- glm(S ~ X + I(X^2), family = binomial(), subset = D == 1)
  }
  
  samp <- predict(sampmod, newdata = data.frame(X = X), type = "response")
  IOW <- (1 - samp)/samp
  
  if (out_misspecify) {
    
    # OLS
    outmod <- glm(Y ~ X, subset = (D == 1 & S == 1), family = gaussian())
    out <- predict(outmod, newdata = data.frame(X = X), type = "response")
    
    # WLS
    outmod.w <- glm(Y ~ X, weights = IOW, subset = (D == 1 & S == 1), family = gaussian)
    out.w <- predict(outmod.w, newdata = data.frame(X = X), type = "response")
    
    # DR Learner approach
    DR <- (Y - out)
    outmod.dr <- glm(DR ~ X, subset = (D == 1 & S == 1), weights = IOW, family = gaussian())
    out.dr <- predict(outmod.dr, newdata = data.frame(X = X), type = "response") + out
    
  } else {
    
    # OLS
    outmod <- glm(Y ~ X + I(X^2), subset = (D == 1 & S == 1), family = gaussian())
    out <- predict(outmod, newdata = data.frame(X = X), type = "response")
    
    # WLS
    outmod.w <- glm(Y ~ X + I(X^2), weights = IOW, subset = (D == 1 & S == 1), family = gaussian)
    out.w <- predict(outmod.w, newdata = data.frame(X = X), type = "response")
    
    # DR Learner approach
    DR <- (Y - out)
    outmod.dr <- glm(DR ~ X + I(X^2), subset = (D == 1 & S == 1), weights = IOW, family = gaussian())
    out.dr <- predict(outmod.dr, newdata = data.frame(X = X), type = "response") + out
    
  }
  
  rmse <- sqrt(mean(Y[D == 0 & S == 0] - out[D == 0 & S == 0])^2)
  rmse.w <- sqrt(mean(Y[D == 0 & S == 0] - out.w[D == 0 & S == 0])^2)
  rmse.dr <- sqrt(mean(Y[D == 0 & S == 0] - out.dr[D == 0 & S == 0])^2)
  
  return(c(rmse, rmse.w, rmse.dr))
  
}

sim_out <- replicate(n.iter, sim_fun(out_misspecify = FALSE, samp_misspecify = TRUE))
rowMeans(sim_out)

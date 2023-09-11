library(mgcv)

n.iter <- 1000

sim_fun <- function(out_misspecify = FALSE, samp_misspecify = FALSE) {
  
  X <- runif(1000, 0, 10)
  V <- runif(1000, 0, 10)
  
  mu <- 1 + X + 0.5*V
  Y <- rnorm(1000, mu, sqrt(X + V))

  eta <- 1.5 - 0.3*X - 0.5*V
  S <- rbinom(1000, 1, plogis(eta))
  D <- rep(c(0,1), each = 500)
  
  if (samp_misspecify == TRUE) {
    sampmod <- glm(S ~ X, family = binomial(), subset = D == 1)
  } else {
    sampmod <- glm(S ~ X + V, family = binomial(), subset = D == 1)
  }
  
  samp <- predict(sampmod, newdata = data.frame(X = X, V = V), type = "response")
  IOW <- (1 - samp)*mean(S)/(samp*(1 - mean(S)))
  
  if (out_misspecify) {
    
    # OLS
    outmod <- glm(Y ~ X, subset = (D == 1 & S == 1), family = gaussian())
    out <- predict(outmod, newdata = data.frame(X = X), type = "response")
    
    # WLS
    outmod.w <- glm(Y ~ X, weights = IOW, subset = (D == 1 & S == 1), family = gaussian)
    out.w <- predict(outmod.w, newdata = data.frame(X = X), type = "response")
    
    # DB Learner
    DB <- (Y - out)
    outmod.db <- glm(DB ~ X, subset = (D == 1 & S == 1), weights = IOW, family = gaussian())
    out.db <- predict(outmod.db, newdata = data.frame(X = X), type = "response") + out
    
  } else {
    
    # OLS
    outmod <- glm(Y ~ X + V, subset = (D == 1 & S == 1), family = gaussian())
    out <- predict(outmod, newdata = data.frame(X = X, V = V), type = "response")
    
    # WLS
    outmod.w <- glm(Y ~ X + V, weights = IOW, subset = (D == 1 & S == 1), family = gaussian)
    out.w <- predict(outmod.w, newdata = data.frame(X = X, V = V), type = "response")
    
    # DB Learner
    DB <- (Y - out)
    outmod.db <- glm(DB ~ X + V, subset = (D == 1 & S == 1), weights = IOW, family = gaussian())
    out.db <- predict(outmod.db, newdata = data.frame(X = X, V = V), type = "response") + out
    
  }
  
  rmse <- sqrt(mean((Y[D == 0 & S == 0] - out[D == 0 & S == 0])^2))
  rmse.w <- sqrt(mean((Y[D == 0 & S == 0] - out.w[D == 0 & S == 0])^2))
  rmse.db <- sqrt(mean((Y[D == 0 & S == 0] - out.db[D == 0 & S == 0])^2))
  
  return(c(rmse, rmse.w, rmse.db))
  
}

sim_out <- replicate(n.iter, sim_fun(out_misspecify = FALSE, samp_misspecify = FALSE))
rowMeans(sim_out)

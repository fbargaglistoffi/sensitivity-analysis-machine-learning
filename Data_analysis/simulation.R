library(mgcv)
library(truncnorm)
library(SuperLearner)

source("~/Github/sensitivity-analysis-machine-learning/Functions/sensitivity_transport.R")

n.iter <- 10

sim_fun <- function(i, out_misspecify = FALSE, samp_misspecify = FALSE, test = seq(0, 0.25, length.out = 21)) {
  
  print(i)
  
  ## simulate data
  X1 <- rnorm(1000, 1, 2)
  X2 <- rnorm(1000, -1, 2)
  X <- data.frame(X1 = X1, X2 = X2)
  
  mu <- 3 + 2*X1 + X2
  Y <- rtruncnorm(1000, a = 0, mean = mu, sd = abs((X1 + X2)))

  eta <- 0.5*X1 - 0.5*X2
  S <- rbinom(1000, 1, plogis(eta))
  
  if (samp_misspecify == TRUE) {
    sampmod <- glm(S ~ X1, family = binomial())
  } else {
    sampmod <- glm(S ~ X1 + X2, family = binomial())
  }
  
  ## fit models
  samp <- predict(sampmod, newdata = data.frame(X = X), type = "response")
  IOW <- (1 - samp)*mean(S)/(samp*(1 - mean(S)))
  
  if (out_misspecify) {
    
    # OLS
    outmod <- glm(Y ~ X1, subset = (S == 1), family = gaussian())
    out <- predict(outmod, newdata = data.frame(X = X), type = "response")
    
    # WLS
    outmod.w <- glm(Y ~ X1, weights = IOW, subset = (S == 1), family = gaussian)
    out.w <- predict(outmod.w, newdata = data.frame(X = X), type = "response")
    
    # DB Learner == WLS
    # DB <- (Y - out)
    # outmod.db <- glm(DB ~ X, subset = (D == 1 & S == 1), weights = IOW, family = gaussian())
    # out.db <- predict(outmod.db, newdata = data.frame(X = X), type = "response") + out
    
  } else {
    
    # OLS
    outmod <- glm(Y ~ X1 + X2, subset = (S == 1), family = gaussian())
    out <- predict(outmod, newdata = data.frame(X = X), type = "response")
    
    # WLS
    outmod.w <- glm(Y ~ X1 + X2, weights = IOW, subset = (S == 1), family = gaussian)
    out.w <- predict(outmod.w, newdata = data.frame(X = X), type = "response")
    
  }
  
  rmse <- sqrt(mean((Y[S == 0] - out[S == 0])^2))
  rmse.w <- sqrt(mean((Y[S == 0] - out.w[S == 0])^2))
  
  sens <- do.call(rbind, lapply(test, dr_sens, S = S, Y = Y, W = X, out = out, samp = samp, q = log, family = gaussian()))
  sens.w <- do.call(rbind, lapply(test, dr_sens, S = S, Y = Y, W = X, out = out.w, samp = samp, q = log, family = gaussian()))
  
  return(list(rmse = rmse, rmse.w = rmse.w, sens = sens, sens.w = sens.w))
  
}

# Run Simulation
sim_out <- lapply(1:n.iter, sim_fun, out_misspecify = FALSE, samp_misspecify = FALSE)

# Report RMSE
mean(do.call(c, lapply(sim_out, function(z) z$rmse)))
mean(do.call(c, lapply(sim_out, function(z) z$rmse.w)))

# Report Sensitivities
eta <- do.call(rbind, lapply(sim_out, function(z) c(z$sens[,1,drop = FALSE]))) # sensitivity parameter
arg <- do.call(rbind, lapply(sim_out, function(z) c(z$sens[,2,drop = FALSE])))
arg.w <- do.call(rbind, lapply(sim_out, function(z) c(z$sens.w[,2,drop = FALSE])))
sens <- cbind(apply(eta, 2, function(x) mean(do.call(c, x), na.rm = T)),
              apply(arg, 2, function(x) mean(do.call(c, x), na.rm = T)), 
              apply(arg.w, 2, function(x) mean(do.call(c, x), na.rm = T)))

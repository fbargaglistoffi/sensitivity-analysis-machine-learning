setwd("~/Dropbox/Projects/FLS-ML")

rm(list=ls())
source("~/Github/sensitivity-analysis-machine-learning/Functions/outliers_bart.R")
source("~/Github/sensitivity-analysis-machine-learning/Functions/sensitivity_bart.R")
source("~/Github/sensitivity-analysis-machine-learning/Functions/calibration.R")
source("~/Github/sensitivity-analysis-machine-learning/Functions/standardize.R")
source("~/Github/sensitivity-analysis-machine-learning/Functions/sensitivity_transport.R")


# data cleaning
library(devtools)
library(tidyverse)
library(haven)
library(dplyr)

# frameworks
library(mice)
library(caret)
library(mlr)
library(osqp)

# ML libraries
library(SuperLearner)
library(ranger)
library(BART)
library(earth)
library(xgboost)

# Plotting
library(ggplot2)
library(knitr)

# clean data
merged_data <- read_csv("Data/PISA_merged_2015_Italy.csv")

pisa_data <- merged_data %>%
  mutate(across(where(is.character), ~ na_if(.x, "99999")),
         across(where(is.numeric), ~ na_if(.x, 99999)))
pisa_data$AGE <- factor(round(pisa_data$AGE))
pisa_data <- pisa_data[which(pisa_data$Region %in% c("Italy: Campania", "Italy: Lombardia")),]
pisa_data$LOMBARDIA <-  as.numeric(pisa_data$Region == "Italy: Lombardia")

vars = c("SC001Q01TA","SC048Q01NA","SC048Q02NA","SC048Q03NA","ST001D01T","ST004D01T","ST013Q01TA",
         "AGE","ISCEDD","ISCEDO","HEDRES","WEALTH","IMMIG","MISCED","FISCED","BMMJ1","BFMJ2",
         "EMOSUPS","REPEAT","OUTHOURS","MMINS","LMINS","PV1MATH","PV1READ","ANXTEST","MOTIVAT",
         "SCHSIZE","CLSIZE","RATCMP1","LEADPD","SCHAUT","EDUSHORT","STRATIO")

imputed <- pisa_data[,vars]
imputed_data <- mice(imputed, method = 'cart', maxit = 5,  m = 1, remove.collinear = FALSE)
demdata <- mice::complete(imputed_data)
datmodel <- cbind(demdata, CNTSTUID = pisa_data$CNTSTUID, 
                  LOMBARDIA = pisa_data$LOMBARDIA, 
                  PV1FLIT = pisa_data$PV1FLIT)
datmodel <- datmodel[complete.cases(datmodel),] # MICE was causing issues so this happened
datmodel[sapply(datmodel, is.character)] <- lapply(datmodel[sapply(datmodel, is.character)], as.factor)

S <- datmodel$LOMBARDIA
n_1 <- sum(S)
X <- model.matrix(~ ., data = datmodel[,vars])[,-1]
X0 <- unname(X[S == 0,])
X1 <- unname(X[S == 1,])

## SuperLearner
samp <- c(SuperLearner(Y = S, X = datmodel[,vars], family = binomial(link = "logit"),
                       SL.library = c("SL.mean","SL.glmnet","SL.ranger"))$SL.predict)

weights <- (1 - samp)*mean(S)/(samp*mean(1 - S))
weights[weights < 0] <- 0
weights[S == 0] <- 1 
cutoff <- quantile(weights[S == 9], 0.99)
weights[weights > cutoff] <- cutoff

## Approximate (or Exact) Quadratic Balancing

# fit_t <- standardize(X = X1, Z = rep(1, nrow(X1)), target = colMeans(X0),
#                      exact_global = FALSE, scale_sample_size = TRUE, lambda = 0)
# weights <- rep(1, nrow(X))
# weights[S == 1] <- n_1*c(fit_t$weights)
# weights[weights < 0] <- 0

## Exact Entropy Balancing

# A <- unname(as.matrix(S*X))
# b <- colMeans(X0)
# fit_t <- cfit(cmat = A, target = n_1*b)
# weights <- c(fit_t$weights)

## Cross-Validation Create 10 equally size folds

# set.seed(42)
# sample <- datmodel[sample(nrow(datmodel)),]
# folds <- cut(seq(1, nrow(sample)), breaks = 10, labels = FALSE)

## No Cross-Validation

sample <- datmodel

# Perform 10 fold cross validation
system.time({
  
  set.seed(42)
  
  # Segment the data by fold using the which() function 
  index <- which(S == 1)
  
  test <- sample[-index, ]
  train <- sample[index, ]
  wts <- weights[index]
  
  x_train <- train[,c(vars)]
  x_test <- test[,c(vars)]
  y_train <- as.vector(train[,c("PV1FLIT")])
  y_test <- as.vector(test[,c("PV1FLIT")])
  
  y <- c(y_train, y_test)
  x <- rbind(x_train, x_test)
  s <- rep(c(1,0), c(nrow(x_train), nrow(x_test)))
  
  # Fit initial model with ranger
  mod_init <- ranger(y ~ ., data = as.data.frame(cbind(y = y_train, x_train)), num.trees = 5000,
                     max.depth = 10, min.node.size = 5, mtry = max(ncol(x_train)/3, 3))
  fit_train_init <- predict(mod_init, data = x_train, type = "response")$predictions
  fit_test_init <- predict(mod_init, data = x_test, type = "response")$predictions
  fit_init <- c(fit_train_init, fit_test_init)
  
  # Fit initial model with xgboost
  # mod_init <- xgboost(label = y_train, dat = model.matrix(~ ., data = x_train), nrounds = 500)
  # fit_train_init <- predict(mod_init, newdata = model.matrix(~ ., data = x_train))
  # fit_test_init <- predict(mod_init, newdata = model.matrix(~ ., x_test))
  
  # Fit weighted model with ranger
  mod_weight <- ranger(y ~ ., data = as.data.frame(cbind(y = y_train, x_train)), num.trees = 10000,
                       case.weights = wts, max.depth = 10, min.node.size = 5, mtry = max(ncol(x_train)/3, 3))
  fit_train_weight <- predict(mod_weight, data = x_train, type = "response")$predictions
  fit_test_weight <- predict(mod_weight, data = x_test, type = "response")$predictions
  fit_weight <- c(fit_train_weight, fit_test_weight)
  
  # Fit weighted model with xgboost
  # mod_weight <- xgboost(label = y_train, dat = model.matrix(~ ., data = x_train), nrounds = 1000, weight = wts)
  # fit_train_weight <- predict(mod_weight, newdata = model.matrix(~ ., data = x_train))
  # fit_test_weight <- predict(mod_weight, newdata = model.matrix(~ ., x_test))
  
  # Fit DR model with ranger
  eps <- (y_train - fit_train_init)
  mod_dr <- ranger(eps ~ ., data = as.data.frame(cbind(eps = eps, x_train)), num.trees = 5000,
                   case.weights = wts, max.depth = 10, min.node.size = 5, mtry = max(ncol(x_train)/3, 3))
  fit_train_dr <- predict(mod_dr, data = x_train, type = "response")$predictions + fit_train_init
  fit_test_dr <- predict(mod_dr, data = x_test, type = "response")$predictions + fit_test_init
  fit_dr <- c(fit_train_dr, fit_test_dr)
  
  # Fit DR model with xgboost
  # eps <- (y_train - fit_train_init)
  # mod_dr <- xgboost(label = eps, dat = model.matrix(~ ., data = x_train), nrounds = 500, weight = wts)
  # fit_train_dr <- predict(mod_dr, newdata = model.matrix(~ ., data = x_train)) + fit_train_init
  # fit_test_dr <- predict(mod_dr, newdata = model.matrix(~ ., x_test)) + fit_test_init
  
  diag_test_weight <- postResample(y_test, fit_test_weight)
  diag_test_init <- postResample(y_test, fit_test_init)  
  diag_test_dr <- postResample(y_test, fit_test_dr)  
  diag_train_weight <- postResample(y_train, fit_train_weight)  
  diag_train_init <- postResample(y_train, fit_train_init)
  diag_train_dr <- postResample(y_train, fit_train_dr)  
  
  # Model Diagnostics.
  RMSE_rf <- c(diag_train_init[1], diag_train_weight[1], diag_train_dr[1],
               diag_test_init[1], diag_test_weight[1], diag_test_dr[1])
  Rsquared_rf <- c(diag_train_init[2], diag_train_weight[2], diag_train_dr[2],
                   diag_test_init[2], diag_test_weight[2], diag_test_dr[2])
  MAE_rf <- c(diag_train_init[3], diag_train_weight[3], diag_train_dr[3],
              diag_test_init[3], diag_test_weight[3], diag_test_dr[3])
  
  # Sensitivity Analysis
  eta <- seq(0,0.25,length.out = 26)
  sens <- do.call(rbind, lapply(eta, dr_sens, S = s, Y = y, W = x, out = fit_init,
                                samp = weights, q = log, family = gaussian()))
  sens.w <- do.call(rbind, lapply(eta, dr_sens, S = s, Y = y, W = x, out = fit_weight,
                                  samp = weights, q = log, family = gaussian()))
  sens.dr <- do.call(rbind, lapply(eta, dr_sens, S = s, Y = y, W = x, out = fit_dr,
                                   samp = weights, q = log, family = gaussian()))
  
  out <- rbind(RMSE_rf, Rsquared_rf, MAE_rf)
  rownames(out) <- c("RMSE", "Rsquared", "MAE")
  colnames(out) <- c("Train - Initial", "Train - Weighted", "Train - DR", 
                     "Test - Initial", "Test - Weighted", "Test - DR")
  
})

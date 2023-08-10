setwd("~/Library/CloudStorage/Dropbox/Projects/FLS-ML")

rm(list=ls())
source("~/Github/sensitivity-analysis-machine-learning/Functions/outliers_bart.R")
source("~/Github/sensitivity-analysis-machine-learning/Functions/sensitivity_bart.R")
source("~/Github/sensitivity-analysis-machine-learning/Functions/calibration.R")
source("~/Github/sensitivity-analysis-machine-learning/Functions/standardize.R")

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

# Approximate (or Exact) Quadratic Balancing

fit_t <- standardize(X = X[S == 1,], Z = S[S == 1], target = colMeans(X[S == 0,]),
                     exact_global = FALSE, scale_sample_size = TRUE, lambda = 0)
weights <- rep(1, nrow(X))
weights[S == 1] <- n_1*c(fit_t$weights)
weights[weights < 0] <- 0

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

## Increase Sample Size

# sample.tmp <- data.frame(datmodel)[S == 1,]
# idx <- sample(1:sum(S), 10000, replace = TRUE)
# sample_s <- sample.tmp[idx,]
# sample <- rbind(data.frame(sample_s), data.frame(datmodel[S == 0,]))
# S <- c(rep(1, times = 10000), rep(0, times = nrow(X0)))
# X <- model.matrix(formula(paste('~', vars, collapse = "+")), data = sample)[,-1]
# 
# fit_t <- standardize(X = X[S == 1,], Z = S[S == 1], target = colMeans(X[S == 0,]),
#                      exact_global = FALSE, scale_sample_size = TRUE, lambda = 1)
# weights <- rep(1, nrow(X))
# weights[S == 1] <- n_1*c(fit_t$weights)
# weights[weights < 0] <- 0

# Perform 10 fold cross validation
system.time({
  
    # Segment the data by fold using the which() function 
    index <- which(S == 1)
    
    test <- sample[-index, ]
    train <- sample[index, ]
    wts <- weights[index]
    
    x_train <- train[,c(vars)]
    x_test <- test[,c(vars)]
    y_train <- as.vector(train[,c("PV1FLIT")])
    y_test <- as.vector(test[,c("PV1FLIT")])

    # Fit Initial Models
    set.seed(42)
    rf_init <- ranger(y ~ ., data = as.data.frame(cbind(y = y_train, x_train)), num.trees = 10000,
                       max.depth = 10, min.node.size = 5, mtry = max(ncol(x_train)/3, 1))

    # Initial Model Predictions
    fit_train_unweight <- predict(rf_init, data = x_train, type = "response")$predictions
    fit_test_unweight <- predict(rf_init, data = x_test, type = "response")$predictions
    # fitted_base <- predict(SuperLearner(Y = y_train, X = x_train,
    #                                     SL.library = c("SL.mean", "SL.glm", "SL.earth",
    #                                                    "SL.xgboost", "SL.ranger")))$pred
    
    # Fit Weighted Models
    rf_new <- ranger(y ~ ., data = as.data.frame(cbind(y = y_train, x_train)), num.trees = 10000,
                      max.depth = 10, min.node.size = 5, mtry = max(ncol(x_train)/3, 1), case.weights = wts)

    # Predict Pseudo Model Fit
    fit_train_weight <- predict(rf_new, data = x_train, type = "response")$predictions
    fit_test_weight <- predict(rf_new, data = x_test, type = "response")$predictions
    
    rf_test_weight <- postResample(y_test, fit_test_weight)
    rf_test_unweight <- postResample(y_test, fit_test_unweight)  
    rf_train_weight <- postResample(y_train, fit_train_weight)  
    rf_train_unweight <- postResample(y_train, fit_train_unweight)
    
    RMSE_rf <- c(rf_train_unweight[1], rf_train_weight[1], rf_test_unweight[1], rf_test_weight[1])
    Rsquared_rf <- c(rf_train_unweight[2], rf_train_weight[2], rf_test_unweight[2], rf_test_weight[2])
    MAE_rf <- c(rf_train_unweight[3], rf_train_weight[3], rf_test_unweight[3], rf_test_weight[3])
    
    out <- rbind(RMSE_rf, Rsquared_rf, MAE_rf)
    rownames(out) <- c("RMSE", "Rsquared", "MAE")
    colnames(out) <- c("Train - Unweighted", "Train - Weighted", "Test - Unweighted", "Test - Weighted")
    
})

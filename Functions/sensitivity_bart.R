#' Function for sensitivity of prediction analysis for BART
#'
#' @param x_train a data frame or matrix of observed predictors in the training set.
#' @param y_train a vector of observed outcomes.
#' @param x_test a data frame or matrix of observed predictors in the test/prediction set.
#' @param nburn the number of iterations used for the burn-in.
#' @param nsamp the number of iterations used to get the results.
#' @param alpha the level of significance for computing differences between posterior predicted values(original model vs augmented model).
#' @return
#' The function returns a list.
#' \code{diff_ppd}: proportion of statistically different posterior predicted values (original model vs augmented model).
#' \code{rmse_original}: estimated RMSE for the original model.
#' \code{rmse_augmented}: estimated RMSE for the augmented model.
#' \code{rsquared_original}: estimated R squared for the original model.
#' \code{rsquared_augmented}: estimated R squared for the augmented model.
#' @export


sensitivity_bart <- function(x_train, y_train, x_test, nburn, nsamp, alpha){
  
  # Upload required packages
  require(BART)
  require(base)
  require(R.utils)
  require(caret)
  require(dplyr)
  require(bayestestR)
  
  # Don't print in-function messages (source: https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html)
  quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  } 
  
  # Initialize objects and parameters
  alpha = as.numeric(alpha)
  seq <- seq(0.1, 0.5, 0.1)
  diff_ppd <- c()
  rmse_original <- c()
  rmse_augmented <- c()
  rsquared_original <- c()
  rsquared_augmented <- c()
  
  # Create progress bar
  pb <- txtProgressBar(min = 0, max = length(seq), style = 3)
  Sys.sleep(0.1)
  
  # Run BART model (Original Model)
  bart <- quiet(wbart(x.train=x_train,
                      y.train=y_train,
                      x.test=x_test,
                      nskip=nburn,
                      ndpost=nsamp))
  
  # Get Draws from the Posterior Distribution
  ppd <- t(apply(bart$yhat.train,
                 1, function(x) rnorm(n = length(x),
                                      mean=x, sd=bart$sigma)))
  ppd_mean <- apply(ppd,2,mean)
  ppd_sd <- apply(ppd,2,sd)
  
  # RMSE
  rmse_original[1] <- postResample(bart$yhat.train.mean, y_train)[1]
  
  # Rsquared
  rsquared_original[1] <- postResample(bart$yhat.train.mean, y_train)[2]
  
  # Get Credibility intervals for Each Unit Level PPD
  conf.int <- as.data.frame(matrix(NA, ncol = 2, nrow = ncol(ppd)))
  for (k in (1:ncol(ppd))){
    conf.int[k, 1:2] <- as.data.frame(ci(ppd[,k], method = "HDI", ci = 0.90))[2:3]
  }
  names(conf.int) <- c("low", "high")
  
  
  for (i in seq){
    
    # Progress bar
    setTxtProgressBar(pb, which(seq==i))
    
    # Build a Correlated New Predictor
    # Generate x1, make sure y is first column, and X is the last
    xy <- cbind(y_train, x1 = rnorm(length(y_train)), data.matrix(x_train))
    
    # Center and scale
    mns <- apply(xy, 2, mean)
    sds <- apply(xy, 2, sd)
    
    xy2 <- sweep(xy, 2, mns, FUN="-")
    xy2 <- sweep(xy2, 2, sds, FUN="/")
    
    # Find existing correlations
    v.obs <- cor(xy2)
    
    # Remove correlation
    xy3 <- xy2 %*% solve(chol(v.obs))
    
    # New correlation
    r <- v.obs
    r[2,] <- r[,2] <- 0
    r[1,2] <- r[2,1] <- i
    diag(r) <- 1
    
    # Use Cholesky decomposition to generate the new matrix
    xy4 <- xy3 %*% chol(r)
    
    # Undo center and scale
    xy4 <- sweep(xy4, 2, sds, FUN="*")
    xy4 <- sweep(xy4, 2, mns, FUN="+")
    
    # Extract x
    x <- as.vector(xy4[,2])
    
    # Run BART model with New Predictor (augmented model)
    x_train_new <- cbind(x_train, x)
    x_test_new <- cbind(x_test, rep(0, nrow(x_test)))
    new_bart <- quiet(wbart(x.train=x_train_new,
                            y.train=y_train,
                            x.test=x_test_new,
                            nskip=nburn,
                            ndpost=nsamp))
    
    # Get Draws from the Posterior Distribution
    
    # For Training Sample
    new_ppd <- t(apply(new_bart$yhat.train,
                       1, function(x) rnorm(n=length(x),
                                            mean = x, sd = new_bart$sigma)))
    new_ppd_mean <- apply(new_ppd,2,mean)
    new_ppd_sd <- apply(new_ppd,2,sd)
    
    # Stat. different predictions
    diff_ppd[which(seq==i)] <- length(which(new_ppd_mean >= conf.int$high | new_ppd_mean <= conf.int$low))/nrow(x_train)
    
    # RMSE
    rmse_augmented[which(seq==i)] <- postResample(new_bart$yhat.train.mean, y_train)[1]
    
    # R squared
    rsquared_augmented[which(seq==i)] <- postResample(new_bart$yhat.train.mean, y_train)[2]
    
    rm(xy, xy2, xy3, xy4, r, mns, sds)
  }
  
  results <- list(diff_ppd, rmse_original, rmse_augmented, rsquared_original, rsquared_augmented)
  names(results) <- c("Proportion of statistically different PPVs", "RMSE original model", "RMSE augmented model", "Rsquared original model", "Rsquared augmented model")
  return(results) 
  close(pb)
}

#' Function for the detection of outliers in the predicted values for BART
#'
#' @param x_train a data frame or matrix of observed predictors in the training set.
#' @param y_train a vector of observed outcomes.
#' @param x_test a data frame or matrix of observed predictors in the test/prediction set.
#' @param nburn the number of iterations used for the burn-in.
#' @param nsamp the number of iterations used to get the results.
#' @param alpha the level of significance.
#' @return
#' The function returns a list.
#' \code{low}: data frame containing the observations with lower-than-the-mean predictions.
#' \code{high}: data frame containing the observations with higher-than-the-mean predictions.
#' @export


outliers_bart <- function(x_train, y_train, x_test, nburn, nsamp, outlier_coeff){
  
  # Upload required package
  require(BART)
  
  # Don't print in-function messages (source: https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html)
  quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  } 
  
  # Initialize parameters
  alpha = as.numeric(outlier_coeff)
  
  # Run BART model
  bart_predictions <- quiet(wbart(x.train=x_train,
                            y.train=y_train,
                            x.test=x_test,
                            nskip=nburn,
                            ndpost=nsamp))
  
  # Get Draws from the Posterior Distribution
  
  # For Test Sample
  ppd_test <- t(apply(bart_predictions$yhat.test,
                      1, function(x) rnorm(n = length(x),
                      mean=x, sd=bart_predictions$sigma)))
  
  # For Training Sample
  ppd_train <- t(apply(bart_predictions$yhat.train,
                       1, function(x) rnorm(n=length(x),
                       mean=x, sd=bart_predictions$sigma)))
  
  # Indicator Variable for Being in the Training Set
  x_test$TRAIN <- rep(1, nrow(x_test))
  x_train$TRAIN <- rep(0, nrow(x_train))
  
  # Get Posterior Mean and SD
  x_train$ppd_mean <- apply(ppd_train,2,mean)
  x_train$ppd_sd <- apply(ppd_train,2,sd)
  
  x_test$ppd_mean <- apply(ppd_test,2,mean)
  x_test$ppd_sd <- apply(ppd_test,2,sd)
  
  y <- y_train
  train <- cbind(y, x_train)
  
  y <- rep(NA, nrow(x_test))
  test <- cbind(y, x_test)
  data <- rbind(train, test)
  colnames(data)[which(names(data) == "y")] <- as.character(as.list(formula)[[2]]) 
    
  # Get Outliers
  
  # Lower Predictions
  low_sd <- data[which(data$ppd_mean <= mean(data$ppd_mean) - outlier_coeff*mean(data$ppd_sd)),]
  low_ad <- data[which(data$ppd_mean <= median(data$ppd_mean) - outlier_coeff*mad(data$ppd_mean)),]
    
  # Higher Predictions
  high_sd <-  data[which(data$ppd_mean >= mean(data$ppd_mean) + outlier_coeff*mean(data$ppd_sd)),]
  high_ad <-  data[which(data$ppd_mean >= median(data$ppd_mean) + outlier_coeff*mad(data$ppd_mean)),]
  
  results <- list(low_sd, low_ad, high_sd, high_ad, data)
  names(results) <- c("Low outliers (sd)", "Low outliers (mad)", "High outliers (sd)", "High outliers (mad)", "Data")
  return(results) 
  
}
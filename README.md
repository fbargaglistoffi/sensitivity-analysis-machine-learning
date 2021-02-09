#  Sensitivity of Machine Learning Predictions

In this repository we provide the code to reproduce the analysis of the paper ["Assessing Sensitivity of Machine Learning Predictions.A Novel Toolbox with an Application to Financial Literacy"](https://arxiv.org/abs/2102.04382) by F.J. Bargagli Stoffi, K. De Beckker, J. E. Maldonado and K. De Witte.

The code is organized as follows:
* <tt>`Data`</tt>: folder containing the .csv data sets used for the analyses;
* <tt>`Data_analysis`</tt>: folder containing the .R, and .Rmd files to reproduce the analyses;
* <tt>`Functions`</tt>: folder containing the .R functions built to implement the analyses.


# Sensitivity function

The <tt>`sensitivity_bart`</tt> function takes as inputs:
* <tt>`x_train`</tt>: a matrix of predictors in the training sample;
* <tt>`y_train`</tt>: a vector of observed values for the outcome in the training sample;
* <tt>`x_test`</tt>: a matrix of predictors in the test sample;
* <tt>`nburn`</tt>: the number of iterations discarded by thealgorithm for the burn-in;
* <tt>`nsamp`</tt>: the number of iterations used by the algorithm  to get the posterior distribution of the outcome;
* <tt>`alpha`</tt>: the level of significance for computing differences between posterior predicted values (original model vs augmented model.

The function returns:
* <tt>`rmse_original`</tt>: the estimated RMSE for the original model;
* <tt>`rmse_augmented`</tt>: the estimated RMSE for the augmented model;
* <tt>`rsquared_original`</tt>: the estimated R2 for the original model,
* <tt>`rsquared_augmented`</tt>: the estimated R2 for the augmented model.


# Example usage

```R
source("sensitivity_bart.R")
sensitivity <- sensitivity_bart(x_train = x_train,
                                y_train = y_train,
                                x_test = x_test,
                                nburn = 1000, 
                                nsamp = 1000,
                                alpha = 0.1)
```



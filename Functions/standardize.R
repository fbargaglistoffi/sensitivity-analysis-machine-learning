# Quadratic Program for Covariate Shift

#' Quadratic program for nonparametric covariate shift
#' @param Y length n vector of observations
#' @param X n x d matrix of covariates
#' @param S n x 1 vector of sample indicators
#' @param kernel kernel function for generating covariance matrix
#' @param lambda regularization hyper parameter, default 0
#' @param eps_abs absolute error tolerance for solver
#' @param eps_rel relative error tolerance for solver
#' @param verbose T/F for whether osqp output should be reported
#' @param ... Extra arguments for osqp solver
#' @export
standardize <- function(Y, X, S, kernel = kernlab::polydot(degree = 2), lambda = 0,
                        eps_abs = 1e-5, eps_rel = 1e-5, verbose = FALSE, ...) {
  
  # ensure that covariate matrices are matrices and get total number of units
  X0.mat <- model.matrix(~ ., data = X[S == 0,])
  X1.mat <- model.matrix(~ ., data = X[S == 1,])
  n <- nrow(X1.mat)
  
  if (is.null(id))
    id <- 1:nrow(X.mat)
  
  kern1 <- kernlab::kernelMatrix(kernel = kernel, x = X1.mat)
  kern0 <- kernlab::kernelMatrix(kernel = kernel, x = X1.mat, y = X0.mat)
  
  # kern1 <- kernlab::kernelMatrix(kernel = kernel, x = W1)
  # kern0 <- kernlab::kernelMatrix(kernel = kernel, x = W1, y = W0)
  
  # construct linear term vector
  q <- -c(rowSums(kern0))
  
  # construct quadratic matrix
  P <- as.matrix(kern1) + diag(1, n) * lambda
  
  # sum to n weights
  A1 <- Matrix::t(rep(1, n))
  l1 <- n
  u1 <- n
  
  # upper and lower bounds of individual weights
  A2 <- Matrix::Diagonal(n)
  l2 <- rep(0, n)
  u2 <- rep(n, n)
  
  A0 <- rbind(A1, A2)
  l0 <- c(l1, l2)
  u0 <- c(u1, u2)
  
  # set optimization settings
  settings <- do.call(osqp::osqpSettings,
                      c(list(verbose = verbose,
                             eps_rel = eps_rel,
                             eps_abs = eps_abs)))
  
  # solve optimization problem
  solution <- osqp::solve_osqp(P = P, q = q, A = A0, l = l0, u = u0, pars = settings)
  weights <- solution$x
  weights[weights < 0] <- .Machine$double.eps
  
  # compute imbalances
  imbalance <- c(colMeans(X0.mat)) - c(t(X1.mat) %*% weights)/n
  names(imbalance) <- colnames(X1.mat)
  
  return(list(weights = weights, imbalance = imbalance, X1 = X1.mat, X0 = X0.mat))
  
}

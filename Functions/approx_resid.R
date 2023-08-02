# Find approximately balancing weights using quadprog
approx.quadprog <- function(cmat, target, zeta = 0.5,
                            allow.negative.weights = FALSE) {
  
  # The system is effectively
  # minimize zeta * delta^2 + (1 - zeta) * ||gamma||^2
  # subject to
  #   sum gamma = 1
  #   delta + (cmat'gamma)_j >= target_j
  #   delta - (cmat'gamma)_j >= -target_j
  #
  # The last two constraints mean that
  # delta = ||cmat'gamma - target||_infty
  
  Dmat = diag(c(zeta, rep(1 - zeta, nrow(cmat))))
  dvec = rep(0, 1 + nrow(cmat))
  Amat = cbind(
    c(0, rep(1, nrow(cmat))),
    rbind(rep(1, ncol(cmat)), cmat ),
    rbind(rep(1, ncol(cmat)), -cmat))
  bvec = c(1, target, -target)
  
  if (!allow.negative.weights) {
    LB = 1/nrow(cmat)/10000
    Amat = cbind(Amat, rbind(rep(0, nrow(cmat)), diag(rep(1, nrow(cmat)))))
    bvec = c(bvec, rep(LB, nrow(cmat)))
  }
  
  balance.soln = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  gamma = balance.soln$solution[-1]
  
  gamma
  
}

svm.bal <- function(X, y, scale = FALSE, cost = 1) {
  
  n <- dim(X)[1]
  tX <- rbind(X, y)
  ty <- c(rep(1, n), -1)
  
  mod <- e1071::svm(x = tX, y = ty, 
                    kernel = "linear", 
                    type = "C-classification", 
                    cost = cost, scale = scale)
  
  coefs <- rep(0,n)
  nsv <- length(mod$index)
  coefs[mod$index[-nsv]] <- mod$coefs[-nsv]
  
  coefs
  
}
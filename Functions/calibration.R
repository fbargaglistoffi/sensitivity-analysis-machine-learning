# Generic Calibration Function

#' Quadratic program for nonparametric covariate shift
#' @param cmat Constraint Matrix
#' @param target n x d matrix of covariates
#' @param base_weights n x 1 vector of initial weights
#' @param coefs_init initial dual variables
#' @param optim_ctrl list of arguments to be passed on to optim
#' @param ... Extra arguments for optim solver
#' @export
#' 
cfit <- function(cmat, target, base_weights = NULL,
                 coefs_init = NULL, optim_ctrl = list(maxit = 500, reltol = 1e-10),
                 ...) {
  
  if (!is.matrix(cmat))
    stop("cmat must be a matrix")
  
  if (!is.vector(target))
    stop("target must be a vector")
  
  if (is.null(base_weights)) {
    base_weights <- rep(1, nrow(cmat))
  } else if (length(base_weights) != nrow(cmat)) { 
    stop("length(base_weights) != sample size")
  }
  
  # initialize coefs
  if (is.null(coefs_init)) {
    coefs_init <- rep(0, times = ncol(cmat)) 
  } else if (length(coefs_init) != ncol(cmat)) {
    stop("length(coefs_init) != ncol(cmat)")
  }
  
  extraArgs <- list(...)
  
  if (length(extraArgs)) {
    
    arg <- names(formals(stats::optim))
    indx <- match(names(extraArgs), arg, nomatch = 0)
    if (any(indx == 0)) 
      stop(paste("Argument", names(extraArgs)[indx == 0], "not matched"))
    
  }
  
  opt <- stats::optim(coefs_init, lagrange_ent, method = "BFGS", cmat = cmat,
                      base_weights = base_weights, target = target,
                      control = optim_ctrl, ...)
  
  converged <- ifelse(opt$convergence == 0, TRUE, FALSE)
  coefs <- opt$par
  weights <- c( base_weights*exp(-cmat %*% coefs) )
  
  out <- list(weights = weights,
              coefs = coefs,
              converged = converged,
              cmat = cmat,
              target = target,
              base_weights = base_weights, 
              optim_ctrl = optim_ctrl)
  
  class(out) <- "cfit"
  return(out)
  
}

lagrange_ent <- function(coefs, cmat, target, base_weights)  {
  
  temp <- sum(base_weights*exp(-cmat %*% coefs))
  out <- temp + sum(target * coefs) 
  
  return(out)
  
}

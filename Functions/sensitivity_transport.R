
# Weird Doubly-Robust Transport
dr_sens <- function(eta, S, Y, W, out, samp, q = I, family = gaussian(), sl.lib = c("SL.mean", "SL.glm")) {
  
  if (length(Y) != length(S))
    stop("length(Y) != length(S)")
  
  if (length(out) != length(S))
    stop("length(out) != length(S)")
  
  if (length(samp) != length(S))
    stop("length(samp) != length(S)")
  
  IOW <- (1 - samp)*mean(S)/(samp*mean(1 - S))
  q <- match.fun(q)
  
  loss <- ((Y - out)/sd(Y - out))^2
  pseudo0 <- exp(eta*q(Y))
  pseudo <- pseudo0*loss
  
  pmod <- SuperLearner(Y = pseudo[S == 1], X = data.frame(W)[S == 1,], family = family, SL.library = sl.lib)
  pmod0 <- SuperLearner(Y = pseudo0[S == 1], X = data.frame(W)[S == 1,], family = family, SL.library = sl.lib)
  
  p <- c(predict(pmod, newdata = data.frame(W))$pred)
  p0 <- c(predict(pmod0, newdata = data.frame(W))$pred)
  theta <- p/p0

  phi <- sum((I(S == 0)*theta + I(S == 1)*(IOW*pseudo0/p0)*(loss - theta)))/sum(S == 0)
  eif <- I(S == 0)*(theta - phi)/mean(S == 0) + I(S == 1)*(IOW*pseudo0/p0)*(loss - theta)
  var <- sum(eif^2)/sum(S == 0)
  
  return(list(eta = eta, loss = phi, loss.var = var))
  
}

base_sens <- function(eta, S, Y, W, out, q = I(), family = family) {
  
  if (length(Y) != length(S))
    stop("length(Y) != length(S)")
  
  if (length(out) != length(S))
    stop("length(out) != length(S)")
  
  q <- match.fun(q)
  
  loss <- (Y - out)^2
  pseudo0 <- exp(eta*q(Y))
  pseudo <- pseudo0*loss
  
  pmod <- SuperLearner(Y = pseudo[S == 1], X = data.frame(W)[S == 1,], family = family, SL.library = sl.lib)
  pmod0 <- SuperLearner(Y = pseudo0[S == 1], X = data.frame(W)[S == 1,], family = family, SL.library = sl.lib)
  
  p <- c(predict(pmod, newdata = data.frame(W))$pred)
  p0 <- c(predict(pmod0, newdata = data.frame(W))$pred)
  theta <- p/p0
  
  phi <- sum(I(S == 0)*theta)/sum(S == 0)
  eif <- (I(S == 0)*(theta - phi)/mean(S == 0))^2
  var <- mean(eif)/length(eif)
  
  return(list(eta = eta, loss = phi, loss.var = var))
  
}


# Weird Doubly-Robust Transport
dr_sens <- function(eta, S, Y, W, out, samp, q = I, family = gaussian()) {
  
  if (length(Y) != length(S))
    stop("length(Y) != length(S)")
  
  if (length(out) != length(S))
    stop("length(out) != length(S)")
  
  if (length(samp) != length(S))
    stop("length(samp) != length(S)")
  
  IOW <- (1 - samp)/samp
  q <- match.fun(q)
  
  loss <- ((Y - out)/sd(Y - out))^2
  pseudo0 <- exp(eta*q(Y))
  pseudo <- pseudo0*loss
  
  pmod0 <- glm(pseudo0 ~ ., data = data.frame(pseudo0 = pseudo0, W), subset = (S == 1), family = family)
  pmod <- glm(pseudo ~ ., data = data.frame(pseudo = pseudo, W), subset = (S == 1), family = family)
  
  p0 <- predict(pmod0, newdata = data.frame(W), type = "response")
  p <- predict(pmod, newdata = data.frame(W), type = "response")
  theta <- p/p0

  phi <- sum((I(S == 0)*theta + I(S == 1)*(IOW*pseudo0/p0)*(loss - theta)))/sum(S == 0)
  eif <- (I(S == 0)*(theta - phi)/mean(S == 0) + I(S == 1)*(IOW*pseudo0/p0)*(loss - theta)/mean(S == 0))^2
  var <- mean(eif)/length(eif)
  
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
  
  pmod0 <- glm(pseudo0 ~ ., data = data.frame(pseudo0 = pseudo0, W), subset = (S == 1), family = family)
  pmod <- glm(pseudo ~ ., data = data.frame(pseudo = pseudo, W), subset = (S == 1), family = family)
  
  p0 <- predict(pmod0, newdata = data.frame(W), type = "response")
  p <- predict(pmod, newdata = data.frame(W), type = "response")
  theta <- p/p0
  
  phi <- sum(I(S == 0)*theta)/sum(S == 0)
  eif <- (I(S == 0)*(theta - phi)/mean(S == 0))^2
  var <- mean(eif)/length(eif)
  
  return(list(eta = eta, loss = phi, loss.var = var))
  
}

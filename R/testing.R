library(mvtnorm)
#' Nonparametric test
#'
#' @param Y
#' @param A
#' @param W
#'
#' @return
#' @export
#'
#' @examples

hteNullTest <- function(Y, A, W) {
  # estimated P(A = 1 | W = w)
  n = length(A)
  prop.reg1 <- do.call(func_1, list(Y=A, X = data.frame(W),
                                    newX = data.frame(W),
                                    family = binomial(),
                                    obsWeights=rep(1,n),
                                    id=1:n))
  pi.hat <- prop.reg1$pred

  # estimated mu
  AW <- cbind(A, data.frame(W))
  if(out.glm) {
    mu.reg  <- glm(Y ~ ., data=AW, family='binomial')
    mu.hat <- mu.reg$fitted.values
    mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1), data.frame(W)), type = 'response')
    mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0), data.frame(W)), type = 'response')
  } else {
    mu.reg <- do.call(func_2, list(Y=Y, X = data.frame(cbind(A, W)),
                                   newX = rbind(data.frame(cbind(A=1, W)), data.frame(cbind(A=0, W))),
                                   family = binomial(),
                                   obsWeights=rep(1,n),
                                   id=1:n))
    mu1.hat <- mu.reg$pred[1:n]
    mu0.hat <- mu.reg$pred[-(1:n)]
    mu.hat <- A * mu1.hat + (1-A) * mu0.hat
  }

  # estimated tau
  tau.hat <- mu1.hat - mu0.hat
  Z.hat <- (2*A - 1) / (A * pi.hat + (1-A) * (1-pi.hat))

  # estimated theta
  gamma.hat <- mean(tau.hat)

  w.ecdf <- function(w) {
    vec.leq <- function(x,y) prod(x <= y)
    return(mean(apply(W, 1, vec.leq, y = w)))
  }
  u.vals <- sort(apply(W, 1, w.ecdf))

  # primitive function
  # n.new * 1 vector
  Gamma.w.vals <- apply(w.vals, 1, function(w0) mean(prod(W <= w0) * tau.hat))
  Omega.w.vals <- Gamma.w.vals - gamma.hat * u.vals

  # nonparametric EIF
  # n * n.new matrix
  eif.Gamma <- apply(w.vals, 1, function(w0) {
    (prod(W <= w0)) * (Z.hat * (Y - mu.hats) + tau.hat)
    - Gamma.w.vals[which(w.vals == w0)]
  })
  eif.Omega <- apply(w.vals, 1, function(w0) {
    (as.numeric(rowSums(W <= w0) == d) - w.ecdf(w0)) * (Z.hat * (Y - mu.hats) + tau.hat - gamma.hat)
    - Omega.w.vals[which(w.vals == w0)]
  })

  # one-step estimators
  # n.new * 1 vector
  Gamma.os.est <- colMeans(eif.Gamma) + Gamma.w.vals
  Omega.os.est <- colMeans(eif.Omega) + Omega.w.vals

  # testing procedure
  # test statistics
  Gamma.stat <- n^1/2*max(abs(Gamma.os.est))
  Omega.stat <- n^1/2*max(abs(Omega.os.est))

  # covariance matrices
  n.new <- length(w.vals)
  Gamma.cov.var <- sapply(1:n.new, function(s) sapply(1:n.new, function(t) {
    mean(eif.Gamma[,s] * eif.Gamma[,t])
  }))
  Omega.cov.var <- sapply(1:length(a.vals), function(s) sapply(1:length(a.vals), function(t) {
    mean(eif.Omega[,s] * eif.Omega[,t])
  }))

  # quantiles
  Gamma.quantile <- qmvnorm(0.95, sigma = Gamma.cov.var, tail = "both")$quantile
  Omega.quantile <- qmvnorm(0.95, sigma = Omega.cov.var, tail = "both")$quantile

  res <- c(Gamma.stat, Omega.stat, Gamma.quantile, Omega.quantile)

  res

}

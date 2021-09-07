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

  # primitive function
  Gamma.hat <- function(w) mean((as.numeric(W <= w)) * tau.hat)
  Omega.hat <- mean()

  # nonparametric EIF
  eif.Gamma <- ((Y - mu.hat) * Z.hat + tau.hat) - Gamma.hat
  eif.Omega <- ((Y - mu.hat) * Z.hat + tau.hat - gamma.hat) - Omega.hat

  # one-step estimators
  Gamma.os.est <- mean(((Y - mu.hat) * Z.hat + tau.hat))
  Omega.os.est <- mean(((Y - mu.hat) * Z.hat + tau.hat - gamma.hat))

  # test statistics
  Gamma.stat <- n^1/2*max(abs(Gamma.os.est))
  Omega.stat <- n^1/2*max(abs(Omega.os.est))

  # covariance matrices
  Gamma.cov.var <- sapply(1:length(a.vals), function(s) sapply(1:length(a.vals), function(t) {
    mean(IF.vals[,s] * IF.vals[,t])
  }))
  Omega.cov.var <- sapply(1:length(a.vals), function(s) sapply(1:length(a.vals), function(t) {
    mean(IF.vals[,s] * IF.vals[,t])
  }))

  # quantiles
  Gamma.quantile <- qmvnorm()
  Omega.quantile <- qmvnorm()

  res <- c(Gamma.stat, Omega.stat, Gamma.quantile, Omega.quantile)

  res

}

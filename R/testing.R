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

hteNullTest <- function(Y, A, W, control = list(), out.glm=FALSE, cov.var=FALSE) {

  # update control parameters
  control <- hte.measure.NullTest.control(control)
  n = length(A)

  if (out.glm){
    prop.reg <- gam(A ~ s(W[,1]) + s(W[,2]) + s(W[,3]), family = 'binomial')
    pi.hat <- prop.reg$fitted.values
    AW <- cbind(A, data.frame(W))
    mu.reg <- glm(Y ~ ., data=AW, family='binomial')
    mu.hat <- mu.reg$fitted.values
    mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1),
                                               data.frame(W)), type = 'response')
    mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0),
                                               data.frame(W)), type = 'response')
    mu.hats <- data.frame(mu1=mu1.hat, mu0=mu0.hat)
    control$pi.hat = pi.hat
    control$mu.hats = mu.hats
  }

  if(control$verbose) cat("Estimating...\n")
  tm0 <- proc.time()
  # estimated propensity
  if(is.null(control$pi.hat)){
    prop.reg <- SuperLearner(Y=A, X = data.frame(W),
                             newX = data.frame(W),
                             SL.library = control$pi.SL.library,
                             family = binomial(),
                             obsWeights=rep(1,n),
                             id=1:n)
    control$pi.hat <- prop.reg$SL.predict
  }

  # estimated outcome regression
  if(is.null(control$mu.hats)){
    AW <- cbind(A, data.frame(W))
    if(length(setdiff(Y, c(0,1))) == 0) {
      family = 'binomial'
    } else {
      family = 'gaussian'
    }
    mu.reg <- SuperLearner(Y=Y, X = data.frame(cbind(A, W)),
                           newX = rbind(data.frame(cbind(A=1, W)), data.frame(cbind(A=0, W))),
                           SL.library = control$mu.SL.library,
                           family = family,
                           obsWeights=rep(1,n),
                           id=1:n)
    control$mu.hats <- data.frame(mu1=mu.reg$SL.predict[1:n], mu0=mu.reg$SL.predict[-(1:n)])
  }

  control$mu.hat <- A * control$mu.hats$mu1 + (1-A) * control$mu.hats$mu0

  # estimated tau
  tau.hat <- control$mu.hats$mu1 - control$mu.hats$mu0
  Z.hat <- (2*A - 1) / (A * control$pi.hat + (1-A) * (1-control$pi.hat))

  # estimated theta
  gamma.hat <- mean(tau.hat)

  w.ecdf <- function(w) {
    vec.leq <- function(x,y) prod(x <= y)
    return(mean(apply(W, 1, vec.leq, y = w)))
  }
  u.vals <- apply(W, 1, w.ecdf)
  # primitive function
  if(control$verbose) cat("Computing Gamma and Omega...\n")
  tm1 <- proc.time()
  # n.new * 1 vector
  w.vals <- W
  Gamma.w.vals <- apply(w.vals, 1, function(w0)
    mean(apply(W, 1, function(x) prod(x <= w0))*tau.hat))
  Omega.w.vals <- Gamma.w.vals - gamma.hat * u.vals

  # nonparametric EIF
  # n * n.new matrix
  vec.eq <- function(x,y) prod(x == y)
  eif.Gamma <- apply(w.vals, 1, function(w0) {
    (apply(W, 1, function(x) prod(x <= w0))) * (Z.hat * (Y - control$mu.hat) + tau.hat) -
      Gamma.w.vals[which(as.logical(apply(W, 1, vec.eq, y = w0)))]
  })
  eif.Omega <- apply(w.vals, 1, function(w0) {
    (apply(W, 1, function(x) prod(x <= w0)) - w.ecdf(w0)) * (Z.hat * (Y - control$mu.hat) + tau.hat - gamma.hat) -
      Omega.w.vals[which(as.logical(apply(W, 1, vec.eq, y = w0)))]
  })

  # one-step estimators
  # n.new * 1 vector
  Gamma.os.est <- colMeans(eif.Gamma) + Gamma.w.vals
  Omega.os.est <- colMeans(eif.Omega) + Omega.w.vals

  # testing procedure
  if(control$verbose) cat("Computing statistics...\n")
  tm2 <- proc.time()
  # test statistics
  Gamma.stat <- n^(1/2)*max(abs(Gamma.os.est))
  Omega.stat <- n^(1/2)*max(abs(Omega.os.est))

  # covariance matrices
  n.new <- dim(w.vals)[1]
  if (cov.var){
    Gamma.cov.var <- sapply(1:n.new, function(s) sapply(1:n.new, function(t) {
      mean(eif.Gamma[,s] * eif.Gamma[,t])
    }))
    Omega.cov.var <- sapply(1:n.new, function(s) sapply(1:n.new, function(t) {
      mean(eif.Omega[,s] * eif.Omega[,t])
    }))

    # quantiles
    tm3 <- proc.time()
    Gamma.epsilon <- rmvnorm(n=control$n.boot, mean=rep(0, n.new), sigma = Gamma.cov.var)
    Gamma.epsilon.stats <- apply(Gamma.epsilon, 1, function(x) {max(abs(x))})
    Omega.epsilon <- rmvnorm(n=control$n.boot, mean=rep(0, n.new), sigma = Omega.cov.var)
    Omega.epsilon.stats <- apply(Omega.epsilon, 1, function(x) {max(abs(x))})
  }else{
    tm3 <- proc.time()
    Gamma.epsilon.stats <- replicate(control$n.boot, max(abs(t(eif.Gamma)%*%rnorm(n.new, 0, 1)/sqrt(n))))
    Omega.epsilon.stats <- replicate(control$n.boot, max(abs(t(eif.Omega)%*%rnorm(n.new, 0, 1)/sqrt(n))))
  }

  Gamma.pvalue <- mean(Gamma.epsilon.stats > Gamma.stat)
  Gamma.quantile <- unname(quantile(Gamma.epsilon.stats, control$conf.level))
  Omega.pvalue <- mean(Omega.epsilon.stats > Omega.stat)
  Omega.quantile <- unname(quantile(Omega.epsilon.stats, control$conf.level))

  ret <- data.frame(type = 'Gamma.stat', stat = Gamma.stat, pvalue = Gamma.pvalue,
                    quantile = Gamma.quantile)
  ret <- rbind(ret,
               data.frame(type = 'Omega.stat', stat = Omega.stat, pvalue = Omega.pvalue,
                          quantile = Omega.quantile))
  tm4 <- proc.time()

  if (control$verbose){
    cat("tm1-tm0 = ", tm1-tm0, "\n")
    cat("tm2-tm1 = ", tm2-tm1, "\n")
    cat("tm3-tm2 = ", tm3-tm2, "\n")
    cat("tm4-tm3 = ", tm4-tm3, "\n")
  }

  ret

}


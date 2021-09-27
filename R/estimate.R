#' Estimator psi(P) and corresponding Confidence Interval for treatment effect over the entire population.
#' Estimator theta(P) and corresponding Confidence Interval for overall measure of the treatment effect heterogeneity.
#'
#' @param A a binary treatment or exposure.
#' @param W a vector of covariates observed prior to A.
#' @param Y the observed outcome.
#' @param control Optional list of control parameters. If \code{control=list()}, default.control.list will be used. See \code{hte.measure.NullTest.control} for details.
#'
#' @return Returns a class of estimator and confidence interval.
#' @export
#'
#' @examples
#' ret1 <- hte.estimator(A, W, Y, control = list(pi.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
#'                                               mu.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
#'                                               conf.int = FALSE, conf.int.type = 'Wald', conf.level = 0.95))
#' ret2 <- hte.estimator(A, W, Y, control = list(pi.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
#'                                               mu.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
#'                                               conf.int = TRUE, conf.int.type = 'Wald', conf.level = 0.95))
#' ret3 <- hte.estimator(A, W, Y, control = list(pi.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
#'                                               mu.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
#'                                               conf.int = TRUE, conf.int.type = 'boot', conf.level = 0.95,
#'                                               n.boot = 500))





htem.estimator <- function(A, W, Y, control = list()){

  control <- hte.measure.NullTest.control(control)
  n = length(A)

  # estimated P(A = 1 | W = w)
  if(is.null(control$pi.hat)){
    prop.reg <- SuperLearner(Y=A, X = data.frame(W),
                              newX = data.frame(W),
                              SL.library = control$pi.SL.library,
                              family = binomial(),
                              obsWeights=rep(1,n),
                              id=1:n)
    control$pi.hat <- prop.reg$SL.predict
  }

  # estimated mu
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

  # 1. estimated psi
  # a) plug-in
  psi.plug.in <- mean(tau.hat^2)
  psi.eif.hat <- 2 * tau.hat * (Y - control$mu.hat) * Z.hat + tau.hat^2 - psi.plug.in
  psi.se <- sd(psi.eif.hat)

  # b) one-step
  psi.one.step.est <- mean(2 * tau.hat * (Y - control$mu.hat) * Z.hat + tau.hat^2)

  # c) truncated
  psi.est <- ifelse(psi.one.step.est>0, psi.one.step.est, psi.plug.in)

  # 2. estimated theta
  gamma.hat <- mean(tau.hat)

  # a) plug-in
  theta.plug.in <- mean((tau.hat-gamma.hat)^2)
  theta.eif.hat <- 2 * (tau.hat-gamma.hat) * (Y - control$mu.hat) * Z.hat + (tau.hat-gamma.hat)^2 - theta.plug.in
  theta.se <- sd(theta.eif.hat)

  # b) one-step
  theta.one.step.est <- mean(2 * (tau.hat-gamma.hat) * (Y - control$mu.hat) * Z.hat + (tau.hat-gamma.hat)^2)

  # c) truncated
  theta.est <- ifelse(theta.one.step.est>0, theta.one.step.est, theta.plug.in)

  ret <- data.frame(type = 'psi.est', est = psi.est, se = psi.se)
  ret <- rbind(ret,
               data.frame(type = 'theta.est', est = theta.est, se = theta.se))

  # confidence interval
  if(control$conf.int){
    if(tolower(control$conf.int.type) == "wald"){
      # Wald-type CI
      ret$ll = ret$est - qnorm(1-(1-control$conf.level)/2) * ret$se / sqrt(n)
      ret$ul = ret$est + qnorm(1-(1-control$conf.level)/2) * ret$se / sqrt(n)
    } else {
      # Bootstrap CI
      boot.ests <- sapply(1:control$n.boot, function(x) {
        boot.inds <- sample(1:n, n, replace=TRUE)

        boot.pi.hat <- control$pi.hat[boot.inds]
        boot.mu.hats <- control$mu.hats[boot.inds,]
        boot.ret <- hte.estimator(A[boot.inds], W[boot.inds], Y[boot.inds],
                             control = list(pi.hat = boot.pi.hat,
                                            mu.hats = boot.mu.hats,
                                            conf.int = FALSE))
        boot.ret$est
      })
      boot.ci <- apply(boot.ests, 1, quantile, c((1-control$conf.level)/2, 1-(1-control$conf.level)/2))
      ret$ll = c(boot.ci[,1][[1]], boot.ci[,2][[1]])
      ret$ul = c(boot.ci[,1][[2]], boot.ci[,2][[2]])
    }
  }
  return(ret)
}

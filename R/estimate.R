#' Estimator psi(P) and corresponding Confidence Interval for treatment effect over the entire population.
#' Estimator theta(P) and corresponding Confidence Interval for overall measure of the treatment effect heterogeneity.
#'
#' @param A a binary treatment or exposure.
#' @param W a vector of covariates observed prior to A.
#' @param Y the observed outcome.
#' @param control Optional list of control parameters. If \code{control=list()}, default control setting will be used.
#' See \code{hte.measure.NullTest.control} for details.
#'
#' @return Returns a class of estimators and confidence intervals.
#' @export
#'
#' @examples
#' ret1 <- htem.estimator(A, W, Y, control = list()) # default control setting-point estimate with only hybrid outputs
#' ret1 <- htem.estimator(A, W, Y, control = list(est.type = list('psi.est'='all', 'theta.est'='hybrid'))) # default control setting-point estimate with additional optional outputs
#' ret1 <- htem.estimator(A, W, Y, control = list(conf.int = TRUE, conf.int.type = 'Wald')) # estimate with wald-type CI
#' ret1 <- htem.estimator(A, W, Y, control = list(conf.int = TRUE, conf.int.type = 'boot', n.boot = 500)) # estimate with Bootstrap CI

htem.estimator <- function(A, W, Y, control = list()){

  # update control parameters
  control <- hte.measure.NullTest.control(control)
  n = length(A)

  # estimated propensity
  if(is.null(control$pi.hat)){
    prop.reg <- SuperLearner(Y=A, X = data.frame(W),
                              newX = data.frame(W),
                              SL.library = control$pi.SL.library,
                              family = binomial(),
                              obsWeights=rep(1,n),
                              id=1:n)
    control$pi.hat <- prop.reg$SL.predict
    if ((min(control$pi.hat)<=0) | (max(control$pi.hat)>=1)) {
      stop("pi.hat <= 0 or pi.hat >= 1")
    }
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

  # 1. estimated psi
  # a) plug-in
  psi.plug.in <- mean(tau.hat^2)
  psi.eif.hat <- 2 * tau.hat * (Y - control$mu.hat) * Z.hat + tau.hat^2 - psi.plug.in
  psi.se <- sd(psi.eif.hat)

  # b) one-step
  psi.one.step.est <- mean(2 * tau.hat * (Y - control$mu.hat) * Z.hat + tau.hat^2)

  # c) hybrid
  psi.est <- ifelse(psi.one.step.est>0, psi.one.step.est, psi.plug.in)

  # 2. estimated theta
  gamma.hat <- mean(tau.hat)

  # a) plug-in
  theta.plug.in <- mean((tau.hat-gamma.hat)^2)
  theta.eif.hat <- 2 * (tau.hat-gamma.hat) * (Y - control$mu.hat) * Z.hat + (tau.hat-gamma.hat)^2 - theta.plug.in
  theta.se <- sd(theta.eif.hat)

  # b) one-step
  theta.one.step.est <- mean(2 * (tau.hat-gamma.hat) * (Y - control$mu.hat) * Z.hat + (tau.hat-gamma.hat)^2)

  # c) hybrid
  theta.est <- ifelse(theta.one.step.est>0, theta.one.step.est, theta.plug.in)

  # estimator
  est.type.names <- names(control$est.type)
  psi.ret <- data.frame()
  theta.ret <- data.frame()
  optional.ret <- data.frame()

  # psi estimator
  if('psi.est' %in% est.type.names){
    psi.ret <- data.frame(type = 'psi.est', est = psi.est, se = psi.se)
    if (control$est.type[["psi.est"]] == 'all'){
      optional.ret <- rbind(optional.ret,
                   data.frame(type = 'psi.plug.in.est', est = psi.plug.in, se = psi.se))
      optional.ret <- rbind(optional.ret,
                   data.frame(type = 'psi.one.step.est', est = psi.one.step.est, se = psi.se))
    }
  }

  # theta estimator
  if('theta.est' %in% est.type.names){
    theta.ret <- data.frame(type = 'theta.est', est = theta.est, se = theta.se)
    if (control$est.type[["theta.est"]] == 'all'){
      optional.ret <- rbind(optional.ret,
                    data.frame(type = 'theta.plug.in.est', est = theta.plug.in, se = theta.se))
      optional.ret <- rbind(optional.ret,
                         data.frame(type = 'theta.one.step.est', est = theta.one.step.est, se = theta.se))
    }
  }
  ret <- rbind(psi.ret, theta.ret)

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
        boot.ret <- htem.estimator(A[boot.inds], W[boot.inds,], Y[boot.inds],
                             control = list(pi.hat = boot.pi.hat,
                                            mu.hats = boot.mu.hats,
                                            conf.int = FALSE))
        boot.ret$ret$est
      })
      boot.ci <- apply(boot.ests, 1, quantile, c((1-control$conf.level)/2, 1-(1-control$conf.level)/2))
      ret$ll = c(boot.ci[,1][[1]], boot.ci[,2][[1]])
      ret$ul = c(boot.ci[,1][[2]], boot.ci[,2][[2]])
    }
  }

  if (length(optional.ret)>0){
    est.ret = list(ret=ret, optional.ret=optional.ret)
  }else{
    est.ret = list(ret=ret)
  }
  return(est.ret)
}

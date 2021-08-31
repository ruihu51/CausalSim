library(plyr)
library(gam)
library(earth)
library(SuperLearner)

#' Estimator psi(P) and corresponding Confidence Interval for treatment effect over the entire population.
#' Estimator theta(P) and corresponding Confidence Interval for overall measure of the treatment effect heterogeneity.
#' Estimator for average treatment effect.
#'
#' @param A a binary treatment or exposure.
#' @param W a vector of covariates observed prior to A.
#' @param Y the observed outcome.
#' @param func_1 superlearner method to estimate P(A = 1 | W = w).
#' @param func_2 superlearner method to estimate mu.
#' @param out.glm If True, estimating mu using glm function.
#'
#' @return Returns a class of estimator and confidence interval.
#' @export
#'
#' @examples ret <- est.psi(A, W, Y, func_1 = "SL.glm", func_2 = "SL.glm")
est.psi <- function(A, W, Y, func_1, func_2, out.glm=FALSE){
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

  # 1. estimated psi
  # a) plug-in
  plug.in <- mean(tau.hat^2)
  eif.hat <- 2 * tau.hat * (Y - mu.hat) * Z.hat + tau.hat^2 - plug.in
  se <- sd(eif.hat)

  ret <- data.frame(type = 'Plug-in', est = plug.in, ll=plug.in - qnorm(.975) * se / sqrt(n), ul=plug.in + qnorm(.975) * se / sqrt(n))

  # b) one-step
  one.step.est <- mean( 2 * tau.hat * (Y - mu.hat) * Z.hat + tau.hat^2)
  ret <- rbind(ret,
               data.frame(type = 'One-step', est = one.step.est, ll=one.step.est - qnorm(.975) * se / sqrt(n), ul=one.step.est + qnorm(.975) * se / sqrt(n)))

  # c) TMLE-new updated
  new.tmle <- ifelse(one.step.est>0, one.step.est, plug.in)
  ret <- rbind(ret,
               data.frame(type = 'TMLE-new', est = new.tmle, ll=new.tmle - qnorm(.975) * se / sqrt(n), ul=new.tmle + qnorm(.975) * se / sqrt(n)))

  # 2. estimated theta
  gamma.hat <- mean(tau.hat)

  # a) plug-in
  plug.in.theta <- mean((tau.hat-gamma.hat)^2)
  se.theta <- sd( 2 * (tau.hat-gamma.hat) * (Y - mu.hat) * Z.hat + (tau.hat-gamma.hat)^2 - plug.in.theta)
  ret <- rbind(ret,
               data.frame(type = 'Plug-in (Theta)', est = plug.in.theta, ll=plug.in.theta - qnorm(.975) * se.theta / sqrt(n), ul=plug.in.theta + qnorm(.975) * se.theta / sqrt(n)))

  # b) one-step
  one.step.est.theta <- mean( 2 * (tau.hat-gamma.hat) * (Y - mu.hat) * Z.hat + (tau.hat-gamma.hat)^2)
  ret <- rbind(ret,
               data.frame(type = 'One-step (Theta)', est = one.step.est.theta, ll=one.step.est.theta - qnorm(.975) * se.theta / sqrt(n), ul=one.step.est.theta + qnorm(.975) * se.theta / sqrt(n)))

  # c) TMLE-new
  new.tmle.theta <- ifelse(one.step.est.theta>0, one.step.est.theta, plug.in.theta)
  ret <- rbind(ret,
               data.frame(type = 'TMLE-new (Theta)', est = new.tmle.theta, ll=new.tmle.theta - qnorm(.975) * se.theta / sqrt(n), ul=new.tmle.theta + qnorm(.975) * se.theta / sqrt(n)))

  # 3. estimated ATE
  # a) plug-in
  plug.in.ate <- mean(tau.hat)
  ate.eif <- (Y - mu.hat) * Z.hat + tau.hat - plug.in.ate
  ret <- rbind(ret, data.frame(type = 'Plug-in (ATE)', est = plug.in.ate, ll = plug.in.ate -1.96 * sd(ate.eif) / sqrt(n),
                               ul = plug.in.ate + 1.96 * sd(ate.eif) / sqrt(n)))

  # b) one-step
  os.ate <- mean((Y - mu.hat) * Z.hat + tau.hat)
  ate.eif <- (Y - mu.hat) * Z.hat + tau.hat - os.ate
  ret <- rbind(ret, data.frame(type = 'One-step (ATE)', est = os.ate, ll = os.ate -1.96 * sd(ate.eif) / sqrt(n),
                               ul = os.ate + 1.96 * sd(ate.eif) / sqrt(n)))

  return(ret)
}



library(ggplot2)

est.psi.plot <- function(ests, plot.type){

  p <- NULL
  est.sum.type = c("Plug-in", "One-step", "TMLE-new")

  if (plot.type == 'density') {
    p <- ggplot(subset(ests, (type %in% est.sum.type))) +
      geom_density(aes(sqrt(n) * (est - psi0), color=type)) +
      facet_wrap(~n, scales='free')
  }

  if (plot.type == 'box') {
    P <- ggplot(subset(ests, (type %in% est.sum.type))) +
      geom_boxplot(aes(as.factor(n), log(1 + sqrt(n) * (est - psi0)), color=type))
  }

  return(p)
}

est.theta.plot <- function(ests, plot.type){

  p <- NULL
  est.sum.type = c("Plug-in (Theta)", "One-step (Theta)", "TMLE-new (Theta)")

  if (plot.type == 'density') {
    p <- ggplot(subset(ests, (type %in% est.sum.type))) +
      geom_density(aes(sqrt(n) * (est - theta0), color=type)) +
      facet_wrap(~n, scales='free')
  }

  if (plot.type == 'box') {
    P <- ggplot(subset(ests, (type %in% est.sum.type))) +
      geom_boxplot(aes(as.factor(n), log(1 + sqrt(n) * (est - theta0)), color=type))
  }

  return(p)
}

est.psi.summary <- function(ests){

  est.sum.type = c("Plug-in", "One-step", "TMLE-new")

  summaries <- ddply(subset(ests, (type %in% est.sum.type)), .(n, type), summarize,
                      na = sum(is.na(est)),
                      coverage = mean(ll <= psi0 & psi0 <= ul, na.rm=TRUE),
                      bias = mean(est - psi0, na.rm=TRUE),
                      var = var(est, na.rm=TRUE),
                      mse = mean((est - psi0)^2, na.rm=TRUE))

  return(summaries)

}

est.theta.summary <- function(ests){

  est.sum.type = c("Plug-in (Theta)", "One-step (Theta)", "TMLE-new (Theta)")

  summaries <- ddply(subset(ests, (type %in% est.sum.type)), .(n, type), summarize,
                     na = sum(is.na(est)),
                     coverage = mean(ll <= theta0 & theta0 <= ul, na.rm=TRUE),
                     bias = mean(est - theta0, na.rm=TRUE),
                     var = var(est, na.rm=TRUE),
                     mse = mean((est - theta0)^2, na.rm=TRUE))

  return(summaries)

}




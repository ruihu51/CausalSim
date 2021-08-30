f <- function(x) log(exp(x) - 1)
finv <- function(x) log(exp(x) + 1)
fprime <- function(x) exp(x) / (exp(x) - 1)

expit <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x) - log(1-x)

get.tmle <- function(Y, A, mu0.hat, mu1.hat, pi.hat, verbose=FALSE, max.iter=30) {

  n <- length(Y)
  mu1.new <- mu1.hat
  mu0.new <- mu0.hat
  mu.new <- A * mu1.new + (1 - A) * mu0.new
  tau.hat <- mu1.hat - mu0.hat
  eif.hat <- 2 * tau.hat * ( (Y - mu1.hat) * A /pi.hat - (Y - mu0.hat) * (1-A) / (1-pi.hat)) + tau.hat^2 - mean(tau.hat^2)
  sigma.hat <- sd(eif.hat)
  g.hat <- A * pi.hat + (1-A) * (1-pi.hat)

  k <- 1
  while(TRUE) {
    if(k > max.iter) {
      break
    }
    mu1.old <- mu1.new
    mu0.old <- mu0.new
    mu.old <- mu.new
    Z.old <- (2*A - 1) * (mu1.old - mu0.old) / g.hat
    D.old <- 2 * (Y - mu.old) * Z.old

    if(abs(mean(D.old)) <= sigma.hat / (log(n) * sqrt(n))) break
    #fluc <- glm(Y ~ 0 + I(A * (mu1.old - mu0.old) / pi.hat) + I((1-A) * (mu1.old - mu0.old) / (1-pi.hat)), offset = logit(mu.old), family='binomial')
    # fluc <- glm(Y ~ 0 + I((2*A - 1) * (mu1.old - mu0.old) / g.hat), offset = logit(mu.old), family='binomial')h
    # eps <- coef(fluc)
    # mu1.new <- expit(logit(mu1.old) + eps * (mu1.old - mu0.old) / pi.hat)
    # mu0.new <- expit(logit(mu0.old) - eps * (mu1.old - mu0.old) / (1-pi.hat))
    #

    eps <- mean((Y - mu.old) * Z.old) / mean(Z.old^2)
    mu1.new <- mu1.old + eps * (mu1.old - mu0.old) / pi.hat
    mu0.new <- mu0.old - eps * (mu1.old - mu0.old) / (1-pi.hat)
    # mu1.new <- expit(logit(mu1.old) + eps[1] * (mu1.old - mu0.old) / pi.hat)
    # mu0.new <- expit(logit(mu0.old) + eps[2] * (mu1.old - mu0.old) / (1-pi.hat))
    mu.new <- A * mu1.new + (1 - A) * mu0.new

    Z.new <- (2*A - 1) * (mu1.new - mu0.new) / g.hat
    # D.new <- 2 * (Y - mu.new) * Z.new
    # nu <- 1 / (pi.hat * (1-pi.hat))
    # mean(D.new)
    # (mean(D.old) / mean(Z.old^2)) * (mean((Y - mu.old)  * Z.old * nu) - .5  *mean(Z.old^2 * nu ) * mean(D.old) / mean(Z.old^2) )
    # (mean(D.old) / mean(Z.old^2)) * (mean((Y - mu.old)  * Z.old * ( nu - mean(Z.old^2/ mean(Z.old^2)  * nu )) ))
    #
    # (mean((Y - mu.old)  * Z.old * ( nu - mean(Z.old^2/ mean(Z.old^2)  * nu )) )) / mean(Z.old^2)
    #
    # (mean((Y - mu.new)  * Z.new * ( nu - mean(Z.new^2/ mean(Z.new^2)  * nu )) )) / mean(Z.new^2)
    #
    # mean(Z.old^2)
    # mean(Z.new^2)
    # mean((mu1.old - mu0.old)^2 * (1 + eps * nu)^2 / g.hat^2)
    #
    # blah <- ( nu - mean(Z.old^2/ mean(Z.old^2)  * nu )) / mean(Z.old^2)
    #
    # mean((Y - mu.old) * Z.old * nu)
    # 4 * mean((Y - mu.old) * Z.old)
    # max(nu) * mean((Y - mu.old) * Z.old)
    # mean()

    k <- k + 1
  }
  if(verbose) {
    if(k <= max.iter) print(paste0("TMLE converged in ", k - 1, " iterations"))
    else print(paste0("TMLE failed to converge in ", max.iter, " iterations"))
  }

  mu1.star <- mu1.new
  mu0.star <- mu0.new
  tau.star <- mu1.star - mu0.star

  tmle <- mean(tau.star^2)
  eif <- 2 * tau.star * ( (Y - mu1.star) * A /pi.hat - (Y - mu0.star) * (1-A) / (1-pi.hat)) + tau.star^2 - tmle
  se <- sd(eif) / sqrt(n)
  # ci <- c(tmle * exp(-qnorm(.975) * se / tmle),tmle * exp(qnorm(.975) * se / tmle))
  #ci <- c(finv(f(tmle) - qnorm(.975) * se * fprime(tmle)), tmle + qnorm(.975) * se)
  #ci <- c(tmle * exp(-qnorm(.975) * se / tmle), tmle + qnorm(.975) * se)
  ci <- tmle + c(-1,1) * qnorm(.975) * se #c(tmle * exp(-qnorm(.975) * se / tmle), tmle + qnorm(.975) * se)

  data.frame(type = 'TMLE', est = tmle, ll=ci[1], ul=ci[2])

}

boot.tmle <- function(Y, A, mu0.hat, mu1.hat, pi.hat, n.boot = 500) {
  n <- length(Y)

  boot.ests <- sapply(1:500, function(n.boot) {
    boot.inds <- sample(1:n, n, replace=TRUE)

    get.tmle(Y = Y[boot.inds], A = A[boot.inds], mu0.hat = mu0.hat[boot.inds], mu1.hat = mu1.hat[boot.inds], pi.hat = pi.hat[boot.inds])['tmle']
  })
  quantile(boot.ests, c(.025, .975))

}

new.tmle <- function(Y, A, mu.hat, mu0.hat, mu1.hat, pi.hat) {
  n <- length(Y)
  tau.hat <- mu1.hat - mu0.hat
  Z.hat <- (2 * A - 1)/ (A * pi.hat + (1-A) * (1-pi.hat))
  inf.fn <- 2 * (Y - mu.hat) * tau.hat * Z.hat + tau.hat^2 - mean(tau.hat^2)
  eps.hat <- ifelse(abs(mean(inf.fn)) < 1 / (sqrt(n) * log(n)),
                    1,
                    mean(Y * tau.hat * Z.hat) / mean(mu.hat * tau.hat * Z.hat))
  mu.star <- mu.hat * eps.hat
  mu1.star <- mu1.hat * eps.hat
  mu0.star <- mu0.hat * eps.hat
  tau.star <- mu1.star - mu0.star
  est <- mean(tau.star^2)
  eif <- 2 * tau.star * (Y - mu.star) * Z.hat + tau.star^2 - est
  se <- sd(eif) / sqrt(n)
  ci <- est + c(-1,1) * qnorm(.975) * se

  data.frame(type = 'TMLE-new', est = est, ll=ci[1], ul=ci[2])

}

W <- matrix(runif(3e5, 0, 1), ncol=3)
pi0 <- function(w) expit(sqrt(abs(w[,1])) + w[,3])
mu0 <- function(a, w) expit(sin(3 * w[,1]) + w[,2]^2 + a * w[,1] + a* sqrt(w[,1]) * w[,2])
#(psi0 <- mean((mu0(1,W) - mu0(0,W))^2))

mu0.null <- function(a, w) expit(sin(3 * w[,1]) + w[,2]^2 + w[,1] + sqrt(w[,1]) * w[,2])
#(psi0 <- mean((mu0.null(1,W) - mu0.null(0,W))^2))

ate0 <- mean(mu0.null(1, W) - mu0.null(0,W))

library(plyr)
library(gam)
library(earth)
library(SuperLearner)

null.sims <- TRUE
null.sims <- FALSE
out.glm <- FALSE


n <- 2000

W <- matrix(runif(n*3, 0, 1), ncol=3)
A <- rbinom(n, size = 1, prob = pi0(W))
if(null.sims) {
  Y <- rbinom(n, size = 1, prob = mu0.null(A, W))
} else {
  Y <- rbinom(n, size = 1, prob = mu0(A, W))
}

est.psi <- function(A, W, Y, func_1, func_2, out.glm=FALSE){
  # estimated P(A = 1 | W = w)
  # prop.reg <- gam(A ~ s(W[,1]) + s(W[,2]) + s(W[,3]), family = 'binomial')
  n = length(A)
  prop.reg1 <- do.call(func_1, list(Y=A, X = data.frame(W),
                             newX = data.frame(W),
                             family = binomial(),
                             obsWeights=rep(1,n),
                             id=1:n))
  # prop.reg1 <- SL.earth(Y=A, X = data.frame(W), newX = data.frame(W), family = binomial(), obsWeights=rep(1,n), id=1:n)
  pi.hat <- prop.reg1$pred
  # pi.hat <- prop.reg$fitted.v

  # estimated mu
  AW <- cbind(A, data.frame(W))
  if(out.glm) {
    mu.reg  <- glm(Y ~ ., data=AW, family='binomial')
    mu.hat <- mu.reg$fitted.values
    mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1), data.frame(W)), type = 'response')
    mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0), data.frame(W)), type = 'response')
  } else {
    # mu.reg <- SL.earth(Y=Y, X = data.frame(cbind(A, W)), newX = rbind(data.frame(cbind(A=1, W)),data.frame(cbind(A=0, W))), family = binomial(), obsWeights=rep(1,n), id=1:n)
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

  # estimated psi
  # plug-in
  plug.in <- mean(tau.hat^2)
  eif.hat <- 2 * tau.hat * (Y - mu.hat) * Z.hat + tau.hat^2 - plug.in
  se <- sd(eif.hat)

  ret <- data.frame(type = 'Plug-in', est = plug.in, ll=plug.in - qnorm(.975) * se / sqrt(n), ul=plug.in + qnorm(.975) * se / sqrt(n))

  # one-step
  one.step.est <- mean( 2 * tau.hat * (Y - mu.hat) * Z.hat + tau.hat^2)
  ret <- rbind(ret,
               data.frame(type = 'One-step', est = one.step.est, ll=one.step.est - qnorm(.975) * se / sqrt(n), ul=one.step.est + qnorm(.975) * se / sqrt(n)))

  # TMLE-new updated
  new.tmle <- ifelse(one.step.est>0, one.step.est, plug.in)
  ret <- rbind(ret,
               data.frame(type = 'TMLE-new', est = new.tmle, ll=new.tmle - qnorm(.975) * se / sqrt(n), ul=new.tmle + qnorm(.975) * se / sqrt(n)))

  # TMLE-new
  # new.tmle <- new.tmle(Y=Y, A=A, mu.hat=mu.hat, mu0.hat=mu0.hat, mu1.hat=mu1.hat, pi.hat=pi.hat)
  # ret <- rbind(ret, new.tmle)

  # estimated theta
  gamma.hat <- mean(tau.hat)

  # plug-in
  plug.in.theta <- mean((tau.hat-gamma.hat)^2)
  se.theta <- sd( 2 * (tau.hat-gamma.hat) * (Y - mu.hat) * Z.hat + (tau.hat-gamma.hat)^2 - plug.in.theta)
  ret <- rbind(ret,
               data.frame(type = 'Plug-in (Theta)', est = plug.in.theta, ll=plug.in.theta - qnorm(.975) * se.theta / sqrt(n), ul=plug.in.theta + qnorm(.975) * se.theta / sqrt(n)))

  # one-step
  one.step.est.theta <- mean( 2 * (tau.hat-gamma.hat) * (Y - mu.hat) * Z.hat + (tau.hat-gamma.hat)^2)
  ret <- rbind(ret,
               data.frame(type = 'One-step (Theta)', est = one.step.est.theta, ll=one.step.est.theta - qnorm(.975) * se.theta / sqrt(n), ul=one.step.est.theta + qnorm(.975) * se.theta / sqrt(n)))

  # TMLE-new
  new.tmle.theta <- ifelse(one.step.est.theta>0, one.step.est.theta, plug.in.theta)
  ret <- rbind(ret,
               data.frame(type = 'TMLE-new (Theta)', est = new.tmle.theta, ll=new.tmle.theta - qnorm(.975) * se.theta / sqrt(n), ul=new.tmle.theta + qnorm(.975) * se.theta / sqrt(n)))

  # estimated ATE
  # plug-in
  plug.in.ate <- mean(tau.hat)
  ate.eif <- (Y - mu.hat) * Z.hat + tau.hat - plug.in.ate
  ret <- rbind(ret, data.frame(type = 'Plug-in (ATE)', est = plug.in.ate, ll = plug.in.ate -1.96 * sd(ate.eif) / sqrt(n),
                               ul = plug.in.ate + 1.96 * sd(ate.eif) / sqrt(n)))

  # one-step
  os.ate <- mean((Y - mu.hat) * Z.hat + tau.hat)
  ate.eif <- (Y - mu.hat) * Z.hat + tau.hat - os.ate
  ret <- rbind(ret, data.frame(type = 'One-step (ATE)', est = os.ate, ll = os.ate -1.96 * sd(ate.eif) / sqrt(n),
                               ul = os.ate + 1.96 * sd(ate.eif) / sqrt(n)))

  return(ret)
}

est.psi.sim <-  function(n_range, j_range, func_1, func_2, null.sims=FALSE){
  ests <- ldply(n_range, function(n) {
    ldply(j_range, function(j) {
      # print(j)
      if(j %% 100 == 0) cat(n, j, '\n')
      seed <- sample(1e3:1e8, 1)
      set.seed(seed)

      W <- matrix(runif(n*3, 0, 1), ncol=3)
      A <- rbinom(n, size = 1, prob = pi0(W))
      if(null.sims) {
        Y <- rbinom(n, size = 1, prob = mu0.null(A, W))
        psi0 <- mean((mu0.null(1,W) - mu0.null(0,W))^2)
        theta0 <- var((mu0.null(1,W) - mu0.null(0,W)))
      } else {
        Y <- rbinom(n, size = 1, prob = mu0(A, W))
        psi0 <- mean((mu0(1,W) - mu0(0,W))^2)
        theta0 <- var((mu0(1,W) - mu0(0,W)))
      }

      ret <- est.psi(A, W, Y, func_1, func_2)

      ret$n <- n
      ret$j <- j
      ret$seed <- seed

      ret$psi0 <- psi0
      ret$theta0 <- theta0

      return(ret)
    })
  })

  return(ests)
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




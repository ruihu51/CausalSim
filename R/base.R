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

# boot.tmle <- function(Y, A, mu0.hat, mu1.hat, pi.hat, n.boot = 500) {
#   n <- length(Y)
#
#   boot.ests <- sapply(1:500, function(n.boot) {
#     boot.inds <- sample(1:n, n, replace=TRUE)
#
#     get.tmle(Y = Y[boot.inds], A = A[boot.inds], mu0.hat = mu0.hat[boot.inds], mu1.hat = mu1.hat[boot.inds], pi.hat = pi.hat[boot.inds])['tmle']
#   })
#   quantile(boot.ests, c(.025, .975))
#
# }

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
# mu0 <- function(a, w) expit(sin(3 * w[,1]) + w[,2]^2 + a * w[,1] + a* sqrt(w[,1]) * w[,2])
# mu0 <- function(a, w) expit(3 * w[,1] + w[,2] + a * w[,1])
mu0 <- function(a, w) expit(3 * w[,1] - 2*w[,2] + 2*a - 4*w[,3])

# mu0.null <- function(a, w) expit(sin(3 * w[,1]) + w[,2]^2 + w[,1] + sqrt(w[,1]) * w[,2])
mu0.null <- function(a, w) expit(3 * w[,1] + w[,2] + w[,1])

ate0 <- mean(mu0.null(1, W) - mu0.null(0,W))


hte.measure.NullTest.control <- function(control){
  control.default = list(pi.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                         mu.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                         conf.int = FALSE, conf.int.type = 'Wald',
                         conf.level = 0.95, n.boot = 500)
  control.names <- names(control)
  if(!is.null(control.names)) {
    for (name in control.names) {
      control.default[[name]] <- control[[name]]
    }
  }
  return(control.default)
}

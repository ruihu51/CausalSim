library(plyr)
library(dplyr)
library(gam)
library(earth)
library(SuperLearner)

#' Function to generate simulated data
#'
#' @return named list containing observed data and the global measures of treatment effect heterogeneity
generate.data <- function(n){
  pi0 <- function(w) expit(sqrt(abs(w[,1])) + w[,3])
  mu0.null <- function(a, w) expit(3 * w[,1] + w[,2] + w[,1])

  W <- matrix(runif(n*3, 0, 1), ncol=3)
  A <- rbinom(n, size = 1, prob = pi0(W))
  Y <- rbinom(n, size = 1, prob = mu0.null(A, W))

  psi0 <- mean((mu0.null(1,W) - mu0.null(0,W))^2)
  theta0 <- var((mu0.null(1,W) - mu0.null(0,W)))

  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}

#' Simulation for treatment effect heterogeneity estimators.
#'
#' @param n_range n range for simulation
#' @param j_range j range for simulation
#' @param func_1 superlearner method to estimate P(A = 1 | W = w).
#' @param func_2 superlearner method to estimate mu.
#' @param null.sims If True, simulate the real data under null hypothesis.
#'
#' @return Returns a list of simulation result.
#' @export
#'
#' @examples
#' ests.sim.1 <- ests.sim(200*c(1:10), 1:500, control = list(conf.int = TRUE), null.sims=FALSE, out.glm=TRUE)
#' ests.sim.2 <- ests.sim(200*c(1:10), 1:500, control = list(conf.int = TRUE, conf.int.type = 'boot'), null.sims=FALSE, out.glm=TRUE)
ests.sim <-  function(n_range, j_range, control, out.glm=TRUE){
  ests <- ldply(n_range, function(n) {
    ldply(j_range, function(j) {
      print(j)
      # if(j %% 100 == 0) cat(n, j, '\n')
      seed <- sample(1e3:1e8, 1)
      set.seed(seed)

      # generate simulated data
      simulated.data <- do.call(generate.data, list(n=n))
      Y <- simulated.data$Y
      A <- simulated.data$A
      W <- simulated.data$W

      # using glm to estimate
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

      # using SuperLearner to estimate
      est.ret <- htem.estimator(A, W, Y, control = control)

      ret <- label.result(est.ret$ret, n, j, seed)
      optional.ret <- NULL
      if(!is.null(est.ret$optional.ret)){
        optional.ret <- label.result(est.ret$optional.ret, n, j, seed)
      }

      parameters.ret <- data.frame(type='pars', n=n, j=j, seed=seed,
                                   psi0=simulated.data$psi0,
                                   theta0=simulated.data$theta0)

      est.sim.ret <- bind_rows(ret, optional.ret, parameters.ret)
      return(est.sim.ret)
    })
  })

  return(ests)
}

testing.sim <-  function(n_range, j_range, control, out.glm=FALSE){
  ret <- data.frame()
  for (n in n_range){
    for (j in j_range){
      # if(j %% 100 == 0)
      cat(n, j, '\n')
      seed <- sample(1e3:1e8, 1)
      set.seed(seed)

      # generate simulated data
      simulated.data <- do.call(generate.data, list(n=n))
      Y <- simulated.data$Y
      A <- simulated.data$A
      W <- simulated.data$W

      if (length(ret)>0){
        ret.tmp <- hteNullTest(Y, A, W, control = control, out.glm = out.glm)

        ret.tmp$n <- rep(n, 2)
        ret.tmp$j <- rep(j, 2)
        ret.tmp$seed <- rep(seed, 2)

        ret.tmp$psi0 <- rep(simulated.data$psi0, 2)
        ret.tmp$theta0 <- rep(simulated.data$theta0, 2)

        ret <- rbind(ret, ret.tmp)
      }
      else{
        ret <- hteNullTest(Y, A, W, control = control, out.glm = out.glm)

        ret$n <- rep(n, 2)
        ret$j <- rep(j, 2)
        ret$seed <- rep(seed, 2)

        ret$psi0 <- rep(simulated.data$psi0, 2)
        ret$theta0 <- rep(simulated.data$theta0, 2)
      }
    }
  }

  return(ret)
}















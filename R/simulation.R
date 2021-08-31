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
#' @examples ests <-  est.psi.sim(200*c(1:2), 1:500, func_1 = "SL.gam", func_2 = "SL.gam", null.sims=FALSE)
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

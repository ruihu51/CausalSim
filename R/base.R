#' Initialize control parameters for hteNullTest
#'
#' This function initializes the control parameters for use in \code{hteNullTest}.
#'
#' @param control Named list containing the control options defined by users.
#' If \code{control=list()}, default control setting will be returned.
#'
#' @return Updated named list containing the control options.
hte.measure.NullTest.control <- function(control){
  control.default = list(pi.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                         mu.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                         est.type = list('psi.est'='hybrid', 'theta.est'='hybrid'),
                         conf.int = FALSE,
                         conf.int.type = 'Wald',
                         conf.level = 0.95,
                         n.boot = 500,
                         verbose = FALSE)
  control.names <- names(control)
  if(!is.null(control.names)) {
    for (name in control.names) {
      control.default[[name]] <- control[[name]]
    }
  }
  return(control.default)
}

# Function to label result
label.result <- function(ret, n, j, seed){
  nrows <- dim(ret)[1]
  ret$n <- rep(n, nrows)
  ret$j <- rep(j, nrows)
  ret$seed <- rep(seed, nrows)
  return(ret)
}

# Function to create repeatable seed
create.seed.dict <- function(x){
  x %>%
    subset(type %in% "Gamma.stat") %>%
    select(c("n", "j", "seed")) %>%
    rename(sample_n=n, j=j, seed=seed)
}

# useful functions
f <- function(x) log(exp(x) - 1)
finv <- function(x) log(exp(x) + 1)
fprime <- function(x) exp(x) / (exp(x) - 1)

expit <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x) - log(1-x)


W <- matrix(runif(3e5, 0, 1), ncol=3)
pi0 <- function(w) expit(sqrt(abs(w[,1])) + w[,3])
# mu0 <- function(a, w) expit(sin(3 * w[,1]) + w[,2]^2 + a * w[,1] + a* sqrt(w[,1]) * w[,2])
# mu0 <- function(a, w) expit(3 * w[,1] + w[,2] + a * w[,1])
mu0 <- function(a, w) expit(3 * w[,1] - 2*w[,2] + 2*a - 4*w[,3])

# mu0.null <- function(a, w) expit(sin(3 * w[,1]) + w[,2]^2 + w[,1] + sqrt(w[,1]) * w[,2])
mu0.null <- function(a, w) expit(3 * w[,1] + w[,2] + w[,1])

ate0 <- mean(mu0.null(1, W) - mu0.null(0,W))

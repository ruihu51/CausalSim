CausalSim
========
The R package `CausalSim` provides simulations for treatment effect heterogeneity estimators.

Installation
---------

### Development version

You can install the **development** version of `CausalSim` package from [GitHub
Repository](https://github.com/ruihu51/CausalSim) with:

```
library(devtools)
devtools::install_github("ruihu51/CausalSim")
```


Usage
-----

### Load the package

``` r
require(CausalSim)
```

### Estimation

``` 
# Generate real data
null.sims <- FALSE

n <- 200
set.seed(86563827)
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

# Estimate
ret1 <- hte.estimator(A, W, Y, control = list(pi.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                                      mu.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                                      conf.int = FALSE, conf.int.type = 'Wald', conf.level = 0.95))

ret2 <- hte.estimator(A, W, Y, control = list(pi.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                                              mu.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                                              conf.int = TRUE, conf.int.type = 'Wald', conf.level = 0.95))

ret3 <- hte.estimator(A, W, Y, control = list(pi.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                                      mu.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                                      conf.int = TRUE, conf.int.type = 'boot', conf.level = 0.95,
                                      n.boot = 500))
```

### Simulation
```
ests <-  est.psi.sim(200*c(1:2), 1:500, func_1 = "SL.gam", func_2 = "SL.gam", null.sims=FALSE)
ests
```

### Summary and Plot
```
est.psi.plot(ests, plot.type='density')
est.theta.plot(ests, plot.type='density')
summaries.psi <- est.psi.summary(ests)
summaries.theta <- est.theta.summary(ests)
```




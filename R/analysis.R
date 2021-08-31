library(ggplot2)

#' Plot of simulation result of estimator for psi
#'
#' @param ests simulation result
#' @param plot.type
#'
#' @return
#' @export
#'
#' @examples est.psi.plot(ests, plot.type='density')
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

#' Plot of simulation result of estimator for theta
#'
#' @param ests simulation result
#' @param plot.type
#'
#' @return
#' @export
#'
#' @examples est.theta.plot(ests, plot.type='density')
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

#' Summary of simulation result of estimator for psi
#'
#' @param ests simulation result
#'
#' @return Returns a data frame of summary of simulation result of estimator for psi
#' @export
#'
#' @examples est.psi.summary(ests)
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

#' Summary of simulation result of estimator for theta
#'
#' @param ests simulation result
#'
#' @return Returns a data frame of summary of simulation result of estimator for theta
#' @export
#'
#' @examples est.theta.summary(ests)
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

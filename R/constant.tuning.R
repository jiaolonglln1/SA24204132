#' @title Tuning lambda assuming GIC to be a constant function in R
#' @description Under the GIC standard, find the best lambda for the lasso regression
#' @param y Response variable which is a n*1 vector
#' @param z Designed n*p matrix of composition as explanatory variables
#' @param mu Penalty parameter avoiding zero-sum constraint from affecting convergence
#' @return A numerical solve making GIC to reach its minimum
#' @examples
#' \dontrun{
#' data(data)
#' y <- data[1]
#' z <- data[2:31]
#' mu <- 1
#' r_clambda(y,z,mu)
#' }
#' @export
r_clambda <- function(y,z,mu){
  y <- as.vector(unlist(y))
  z <- as.matrix(z)
  n <- nrow(z); p <- ncol(z)
  beta_hat <- function(lambda){r_cdmlasso(y,z,lambda,mu)}
  sigma2_hat <- function(lambda){(norm((y-(z %*% beta_hat(lambda))),type = '2')^2)/n}
  s <- function(lambda){sum(beta_hat(lambda) != 0)}
  gic <- function(lambda){log(sigma2_hat(lambda))+(s(lambda)-1)*(log(log(n))/n)*log(max(n,p))}
  lambda_better <- optimize(gic,interval = c(-5,5),maximum = F)$minimum
  return(lambda_better)
}

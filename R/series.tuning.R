#' @title Tuning lambda using a prepared series within GIC standard in R
#' @description Under the GIC standard, find the best lambda for the lasso regression
#' @param y Response variable which is a n*1 vector
#' @param z Designed n*p matrix of composition as explanatory variables
#' @param mu Penalty parameter avoiding zero-sum constraint from affecting convergence
#' @param low Lower bound of the interval of prepared lambda
#' @param up Upper bound of the interval of prepared lambda
#' @param by Step size of the interval of prepared lambda
#' @return The best point within the series to minimum GIC
#' @examples
#' \dontrun{
#' data(data)
#' y <- data[1]
#' z <- data[2:31]
#' mu <- 1
#' r_slambda(y,z,mu,low = 0,up = 2)
#' }
#' @export
r_slambda <- function(y,z,mu,low = 0, up = 5, by = 0.1){
  y <- as.vector(unlist(y))
  z <- as.matrix(z)
  n <- nrow(z); p <- ncol(z)
  lambdapool <- seq(low,up,by = by)
  gic <- numeric(length(lambdapool))
  for (i in 1:length(lambdapool)) {
    beta_hat <- r_cdmlasso(y,z,lambdapool[i],mu)
    sigma2_hat <- (norm((y-(z %*% beta_hat)),type = '2')^2)/n
    s <- sum(beta_hat != 0)
    gic[i] <- log(sigma2_hat)+(s-1)*(log(log(n))/n)*log(max(n,p))
  }
  lambda <- lambdapool[which(gic == min(gic))[1]]
  return(lambda)
}


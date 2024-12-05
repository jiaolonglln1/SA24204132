#' @title A Coordinate Descend Algorithm for Compositional Data Lasso problem in R
#' @description This function proposed a coordinate descend algorithm in high dimensional lasso regression problems.
#' @param y Response variable which is a n*1 vector
#' @param z Designed n*p matrix of composition as explanatory variables
#' @param lambda Regularization parameter in Lasso problem
#' @param mu Penalty parameter avoiding zero-sum constraint from affecting convergence(suggested to be 1)
#' @param tol Iteration stop condition: Tolerance
#' @param max_iter Iteration stop condition: Max iterations
#' @return A p*1 vector of regression coefficient estimator
#' @examples
#' \dontrun{
#' data(data)
#' y <- data[1]
#' z <- data[2:31]
#' mu <- 1; lambda <- 0.5
#' r_cdmlasso(y,z,lambda,mu)
#' }
#' @export
r_cdmlasso <- function(y,z,lambda,mu,tol=1e-5,max_iter=1000){
  y <- as.vector(unlist(y))
  z <- as.matrix(z)
  #定义算法中各变量
  n <- nrow(z); p <- ncol(z)
  beta <- rep(0,p); alpha <- 0
  #v <- (unlist(lapply(1:p, function(k) norm(z[,k],type = '2')))^2)/n
  v <- numeric(p)
  for (i in 1:p) {
    v[i] <- (norm(z[,i],type = '2')^2)/n
  }

  #迭代算法
  iter_out <- 0
  while (1) {
    beta_out <- beta
    alpha_out <- alpha
    iter_in <- 0
    while (1) {
      beta_in <- beta_out
      for (j in 1:p) {
        t <- (1/n)*z[,j] %*% (y-z[,-j] %*% beta_in[-j])-mu*(sum(beta_in[-j])+alpha)
        beta_in[j] <- (1/(v[j]+mu))*soft_thr(t,lambda)
        iter_in <- iter_in+1
      }
      if (max(abs(beta_in-beta_out))<=tol | iter_in>max_iter) {
        beta_out <- beta_in
        break
      }
    }
    alpha_out <- alpha_out+sum(beta_out)
    iter_out <- iter_out+1
    if (max(abs(beta_out-beta))<tol | iter_out>max_iter) {
      return(beta_out)
    }else{
      beta <- beta_out
      alpha <- alpha_out
    }
  }
}

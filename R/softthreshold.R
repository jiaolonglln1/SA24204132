#' @title Soft threshold function
#' @param t A specific form in the operation term
#' @param lambda Regularization parameter in Lasso problem
#' @return The Lasso Soft-threshold operation term computed by samples
#' @export
soft_thr <- function(t,lambda){
  sign(t)*max(abs(t)-lambda,0)
}

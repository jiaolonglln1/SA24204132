#' @title A Treated Compositional Dataset For Demonstration
#' @docType data
#' @name data
#' @description A dataset used to benchmark \code{r_cdmlasso} and \code{c_cdmlasso}
#' @examples
#' \dontrun{
#' data(data)
#' y <- as.vector(unlist(data[1]))
#' z <- as.matrix(data[2:31])
#' mu <- 1; lambda <- 0.5
#' microbenchmark::microbenchmark(R=r_cdmlasso(y,z,lambda, mu),C=c_cdmlasso(y,z,lambda, mu))
#' }
NULL

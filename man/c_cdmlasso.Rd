% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{c_cdmlasso}
\alias{c_cdmlasso}
\title{A Coordinate Descend Algorithm for Compositional Data Lasso problem in Rcpp}
\usage{
c_cdmlasso(y, z, lambda, mu, tol = 1e-04, max_iter = 1000L)
}
\arguments{
\item{y}{Response variable which is a n*1 vector}

\item{z}{Designed n*p matrix of composition as explanatory variables}

\item{lambda}{Regularization parameter in Lasso problem}

\item{mu}{Penalty parameter avoiding zero-sum constraint from affecting convergence(suggested to be 1)}

\item{tol}{Iteration stop condition: Tolerance}

\item{max_iter}{Iteration stop condition: Max iterations}
}
\description{
Specific information can be found in 'r_cdmlasso' function(Pay attention that 'c_cdmlasso' only receive numeric vector and matrix as input which is different from 'r_cdmlasso')
}
\examples{
\dontrun{
data(data)
y <- as.vector(unlist(data[1]))
z <- as.matrix(data[2:31])
mu <- 1; lambda <- 0.5
c_cdmlasso(y,z,lambda,mu)
}
}

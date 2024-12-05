#include <Rcpp.h>
using namespace Rcpp;

double c_soft_thr(double t, double lambda) {
    return (t > lambda) ? (t - lambda) : ((t < -lambda) ? (t + lambda) : 0.0);
}

NumericVector m_mult_v(NumericMatrix A, NumericVector b) {
    int nrow = A.nrow();
    int ncol = A.ncol();
    if (ncol != b.size()) {
        stop("Number of columns in matrix must equal the length of the vector.");
    }
    NumericVector result(nrow);
    for (int i = 0; i < nrow; i++) {
        result[i] = 0;
        for (int j = 0; j < ncol; j++) {
            result[i] += A(i, j) * b[j];
        }
    }
    return result;
}

double v_mult_v(NumericVector row, NumericVector col) {
    if (row.size() != col.size()) {
        stop("The dimensions of the row and column vectors must be the same.");
    }
    double result = 0;
    for (int i = 0; i < row.size(); i++) {
        result += row[i] * col[i];
    }
    return result;
}

NumericMatrix remove_column(NumericMatrix mat, int col_to_remove) {
    int nrow = mat.nrow();
    int ncol = mat.ncol();
    NumericMatrix result(nrow, ncol - 1);
    int new_col = 0;
    for (int j = 0; j < ncol; j++) {
        if (j != col_to_remove) {
            result(_, new_col) = mat(_, j);
            new_col++;
        }
    }
    return result;
}

NumericVector remove_row(NumericVector vec, int row_to_remove) {
    int n = vec.size();
    NumericVector result(n - 1);
    int new_index = 0;
    for (int i = 0; i < n; i++) {
        if (i != row_to_remove) {
            result[new_index] = vec[i];
            new_index++;
        }
    }
    return result;
}

//' @title A Coordinate Descend Algorithm for Compositional Data Lasso problem in Rcpp
//' @description Specific information can be found in 'r_cdmlasso' function(Pay attention that 'c_cdmlasso' only receive numeric vector and matrix as input which is different from 'r_cdmlasso')
//' @param y Response variable which is a n*1 vector
//' @param z Designed n*p matrix of composition as explanatory variables
//' @param lambda Regularization parameter in Lasso problem
//' @param mu Penalty parameter avoiding zero-sum constraint from affecting convergence(suggested to be 1)
//' @param tol Iteration stop condition: Tolerance
//' @param max_iter Iteration stop condition: Max iterations
//' @examples
//' \dontrun{
//' data(data)
//' y <- as.vector(unlist(data[1]))
//' z <- as.matrix(data[2:31])
//' mu <- 1; lambda <- 0.5
//' c_cdmlasso(y,z,lambda,mu)
//' }
//' @export
// [[Rcpp::export]]
NumericVector c_cdmlasso(NumericVector y, NumericMatrix z, double lambda, double mu, double tol = 1e-4, int max_iter = 1000) {
    int n = y.size();
    int p = z.ncol();
    NumericVector beta(p, 0.0);
    double alpha = 0.0;
    NumericVector v(p);
    for (int i = 0; i < p; i++) {
        v[i] = sum(pow(z(_, i), 2)) / n;
    }
    for (int iter_out = 0; iter_out < max_iter; iter_out++) {
        NumericVector beta_out = beta;
        double alpha_out = alpha;
        for (int iter_in = 0; iter_in < max_iter; iter_in++) {
            NumericVector beta_in = beta_out;
            for (int j = 0; j < p; j++) {
                NumericVector z_j = z(_, j);
                NumericMatrix z_minus_j = remove_column(z, j);
                NumericVector beta_minus_j = remove_row(beta_in, j);
                double t = (1.0 / n) * v_mult_v(z_j, (y - m_mult_v(z_minus_j, beta_minus_j))) - mu * (sum(beta_minus_j) + alpha);
                beta_in[j] = (1.0 / (v[j] + mu)) * c_soft_thr(t, lambda);
            }
            if (max(abs(beta_in - beta_out)) <= tol) {
                beta_out = beta_in;
                break;
            }
        }
        alpha_out += sum(beta_out);
        if (max(abs(beta_out - beta)) < tol) {
            return beta_out;
        }
        else {
            beta = beta_out;
            alpha = alpha_out;
        }
    }
    return beta;
}

double sigma2_hat(NumericVector y, NumericMatrix z, NumericVector beta_hat) {
    int n = y.size();
    NumericVector residual = y - m_mult_v(z, beta_hat);
    return sum(pow(residual, 2)) / n;
}

int s(NumericVector beta_hat) {
    int non_zero_count = 0;
    for (int i = 0; i < beta_hat.size(); i++) {
        if (beta_hat[i] != 0) {
            non_zero_count++;
        }
    }
    return non_zero_count;
}

double gic(NumericVector y, NumericMatrix z, double lambda, double mu, int n, int p) {
    NumericVector beta_hat = c_cdmlasso(y, z, lambda, mu);
    double sigma2 = sigma2_hat(y, z, beta_hat);
    int s_lambda = s(beta_hat);
    return log(sigma2) + (s_lambda - 1) * (log(log(n)) / n) * log(std::max(n, p));
}

double gic_grad(NumericVector y, NumericMatrix z, double lambda, double mu, int n, int p) {
    double epsilon = 1e-6;
    double gic_plus = gic(y, z, lambda + epsilon, mu, n, p);
    double gic_minus = gic(y, z, lambda - epsilon, mu, n, p);
    return (gic_plus - gic_minus) / (2 * epsilon);  // Numerical gradient using central difference
}

//' @title Tuning lambda assuming GIC to be a constant function in Rcpp
//' @description Specific information can be found in 'r_clambda' function(It has the same problem as 'c_cmdlasso' which means you need to make sure that y and z is numeric before using)
//' @param y Response variable which is a n*1 vector
//' @param z Designed n*p matrix of composition as explanatory variables
//' @param mu Penalty parameter avoiding zero-sum constraint from affecting convergence(suggested to be 1)
//' @param learning_rate Learning rate for upgrading grad function
//' @param max_iter Iteration stop condition: Max iterations
//' @param tol Iteration stop condition: Tolerance
//' @examples
//' \dontrun{
//' data(data)
//' y <- as.vector(unlist(data[1]))
//' z <- as.matrix(data[2:31])
//' mu <- 1
//' c_clambda(y,z,mu)
//' }
//' @export
// [[Rcpp::export]]
double c_clambda(NumericVector y, NumericMatrix z, double mu, double learning_rate = 0.01, int max_iter = 1000, double tol = 1e-6) {
    int n = y.size();
    int p = z.ncol();
    double lambda = 0.8;
    for (int iter = 0; iter < max_iter; iter++) {
        double grad = gic_grad(y, z, lambda, mu, n, p);
        lambda -= learning_rate * grad;
        if (std::abs(grad) < tol) {
            break;
        }
    }
    return lambda;
}

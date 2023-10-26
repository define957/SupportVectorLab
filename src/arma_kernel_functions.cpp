#include <RcppArmadillo.h>

using namespace std;
using namespace arma;

//[[Rcpp::export]]
SEXP cpp_rbf_kernel(arma::mat& x1, arma::mat& x2, double gamma, bool symmetric = 1) {
  int n1 = x1.n_rows;
  arma::mat e1;
  arma::mat norms1 = sum(square(x1), 1) ;
  arma::mat norms2;
  arma::mat e2; e2.ones(1, n1);
  if (symmetric == FALSE) {
    int n2 = x2.n_rows;
    norms2 = sum(square(x2), 1) ;
    e1.ones(1, n2);
  } else {
    norms2 = norms1;
    e1 = e2;
  }
  return Rcpp::wrap(exp(gamma*(-norms1*e1 - e2.t()*norms2.t() + 2*(x1 * x2.t()))));
}

//[[Rcpp::export]]
SEXP cpp_linear_kernel(arma::mat& x1, arma::mat& x2) {
  return Rcpp::wrap(x1 * x2.t());
}

//[[Rcpp::export]]
SEXP cpp_poly_kernel(arma::mat& x1, arma::mat& x2, double degree, double coef0) {
  return Rcpp::wrap(arma::pow(x1 * x2.t() + coef0, degree));
}

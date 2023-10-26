#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
SEXP cpp_chol_solve (arma::mat A, arma::mat b) {
  arma::mat res = arma::solve(A, b, arma::solve_opts::likely_sympd);
  return Rcpp::wrap(res);
}

// [[Rcpp::export]]
SEXP cpp_inverse_spd(arma::mat A) {
  arma::mat res = arma::inv_sympd(A);
  return Rcpp::wrap(res);
}

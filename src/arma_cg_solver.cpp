#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec cpp_conjugate_gradient_method(arma::mat A, arma::vec b, arma::vec x,
                                        const double eps, const unsigned int max_steps) {

  arma::vec rk = b - A * x;
  arma::vec pk = rk;
  double rk2, alphak, betak;

  unsigned int t;

  for (t = 1; t <= max_steps; t++) {
    rk2 = arma::dot(rk, rk);
    arma::vec Apk = A * pk;
    alphak = rk2 / arma::dot(pk, Apk);

    x += alphak * pk;
    arma::vec rk_1 = rk - alphak * Apk;

    if (arma::norm(rk_1, 2) < eps) {
      Rcpp::Rcout << "converge after " << t << " steps.\n";
      break;
    } else {
      betak = arma::dot(rk_1, rk_1) / rk2;
      pk = rk_1 + betak * pk;
    }
    rk = rk_1;
  }

  if (t > max_steps) {
    Rcpp::warning("Conjugate gradient did not converge within %d steps.", max_steps);
  }

  return x;
}

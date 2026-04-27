#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List arma_lsvm_optimizer(const arma::mat& H_left, const arma::mat& H_right,
                               const double a,
                               const arma::mat& A, const arma::mat& B,
                               const arma::vec& q,
                               double eps, int max_steps, arma::vec u,
                               double beta) {
  arma::vec Hu = H_left * (H_right * u) + a*u;
  arma::vec zero_vec = arma::zeros<arma::vec>(q.n_elem);
  arma::vec u_new(q.n_elem);

  int t;
  for (t = 0; t < max_steps; ++t) {
    arma::vec z = (q + arma::max(Hu - q - beta * u, zero_vec));
    arma::vec v = B * z;
    u_new       = a * z - A * v;
    if (arma::norm(u_new - u, 2) < eps) {
      u  = u_new;
      Hu = H_left * (H_right * u) + a*u;
      break;
    }
    u  = u_new;
    Hu = H_left * (H_right * u) + a*u;
  }

  double obj_val = 0.5 * arma::as_scalar(u.t() * Hu) - arma::as_scalar(q.t() * u);
  Rcpp::List res = Rcpp::List::create(Rcpp::Named("x") = u,
                                      Rcpp::Named("iterations") = t,
                                      Rcpp::Named("objective.value") = obj_val);
  return res;
}

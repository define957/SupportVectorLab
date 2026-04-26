#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List arma_lsvm_optimizer(const arma::mat& H, const arma::mat& Hinv,
                               const arma::vec& q,
                               double eps, int max_steps, arma::vec u,
                               double beta) {
  arma::vec Hu = H * u;
  arma::vec Hinv_q = Hinv * q;
  arma::vec zero_vec = arma::zeros<arma::vec>(q.n_elem);
  arma::vec u_new(q.n_elem);

  int t;
  for (t = 0; t < max_steps; ++t) {
    u_new = Hinv_q + Hinv * arma::max(Hu - q - beta * u, zero_vec);
    if (arma::norm(u_new - u, 2) < eps) {
      u  = u_new;
      Hu = H * u;
      break;
    }
    u  = u_new;
    Hu = H * u;
  }

  double obj_val = 0.5 * arma::as_scalar(u.t() * Hu) - arma::as_scalar(q.t() * u);
  Rcpp::List res = Rcpp::List::create(Rcpp::Named("x") = u,
                                      Rcpp::Named("iterations") = t,
                                      Rcpp::Named("objective.value") = obj_val);
  return res;
}

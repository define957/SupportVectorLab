#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
Rcpp::List cpp_clip_dcd_optimizer(arma::mat H, arma::mat q,
                                  arma::mat lb, arma::mat ub,
                                  double eps, unsigned int max_steps,
                                  arma::mat u){
  unsigned int n = H.n_rows;

  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int max_idx = 0;
  double L_val_max;
  double lambda_max;
  double lambda_opt;
  arma::mat numerator(n, 1);
  arma::vec L_idx_val(n);
  arma::vec L_val(n);
  arma::mat Hu(n, n);
  Hu = H * u;
  arma::vec Huik(n);
  arma::mat Hui(n, n);
  arma::uvec idx;
  arma::uvec unique_idx;
  arma::mat diagH = H.diag();
  arma::mat ut = u.t();
  for (i = 0; i < n; i++) {
    Hui(i, span::all) = H.row(i) % ut;
  }
  for(i = 0; i < max_steps; i++){
    numerator = q - Hu;
    L_idx_val = numerator / diagH;
    L_val = numerator % L_idx_val;
    unique_idx = find((u > lb && L_idx_val < 0) || (u < ub && L_idx_val > 0));
    if(unique_idx.n_elem == 0){
      break;
    }
    L_val_max = as_scalar(max(L_val(unique_idx)));
    if(L_val_max < eps){
      break;
    }
    idx = find(L_val == L_val_max);
    for (j = 0; j < idx.n_elem; j++) {
      max_idx = as_scalar(idx(j));
      lambda_max = as_scalar(L_idx_val(max_idx));
      lambda_opt = max(as_scalar(lb(max_idx) - u(max_idx)),
                       min(lambda_max, as_scalar(ub(max_idx) - u(max_idx))));
      if (lambda_opt != 0) {
        break;
      }
    }
    u(max_idx) = u(max_idx) + lambda_opt;
    Huik = H.col(max_idx)*u(max_idx);
    Hu = Hu - Hui.col(max_idx) + Huik;
    Hui(span::all, max_idx) = Huik;
  }
  double obj_val = as_scalar(0.5 * u.t() * H * u - q.t() * u);

  List res = List::create(Named("x") = u,
                          Named("iterations") = i+1,
                          Named("objectiv.value") = obj_val);
  return res;
}

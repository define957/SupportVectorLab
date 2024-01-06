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
  unsigned int max_idx_temp = 0;
  double L_val_max;
  // double lambda_max;
  double lambda_opt = 0;
  double lambda_opt_temp = 0;
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
  arma::vec lambda_max_list;
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
    lambda_max_list = L_idx_val(idx);
    lambda_opt = 0;
    max_idx = 0;
    for (j = 0; j < idx.n_elem; j++) {
      max_idx_temp = as_scalar(idx(j));
      lambda_opt_temp = max(as_scalar(lb(max_idx_temp) - u(max_idx_temp)),
                        min(lambda_max_list(j),
                            as_scalar(ub(max_idx_temp) - u(max_idx_temp))));
      if (abs(lambda_opt) < abs(lambda_opt_temp)) {
        max_idx = max_idx_temp;
        lambda_opt = lambda_opt_temp;
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

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
Rcpp::List cpp_clip_dcd_optimizer(const arma::mat& H, const arma::vec& q,
                                  const arma::vec& lb, const arma::vec& ub,
                                  const double eps, unsigned int max_steps,
                                  arma::vec u) {

  const unsigned int n                    = H.n_rows;
  const arma::vec    inv_diagH            = 1.0 / H.diag();
  unsigned int       iter                 = 0;
  arma::vec          numerator            = q - H * u;
  arma::vec          ub_u                 = ub - u;
  arma::vec          lb_u                 = lb - u;
  arma::uvec         feasible_indices(n);

  for(iter = 0; iter < max_steps; ++iter) {

    unsigned int count = 0;
    for(unsigned int i = 0; i < n; ++i) {

      double nu   = numerator(i);
      double lb_i = lb_u(i);
      double ub_i = ub_u(i);

      if(std::fabs(nu) < 1e-12) {
        continue;
      }
      if((lb_i < 0 && nu < 0) ||
         (ub_i > 0 && nu > 0)) {
        feasible_indices(count++) = i;
      }
    }

    if(count == 0){break;}

    double       max_decrease = 0.0;
    unsigned int max_count    = 0;
    unsigned int best_idx     = 0;

    for(unsigned int k = 0; k < count; ++k) {
      unsigned int idx      = feasible_indices(k);
      double       nu       = numerator(idx);
      double       decrease = (nu * nu) * inv_diagH(idx);

      if(decrease > max_decrease) {
        max_decrease        = decrease;
        max_count           = 1;
        feasible_indices(0) = idx;
        best_idx            = idx;
      } else if(decrease == max_decrease) {
        feasible_indices(max_count++) = idx;
      }
    }

    if(max_decrease < eps){break;}

    if(max_count > 1) {
      double max_step_abs = 0.0;
      for(unsigned int k = 0; k < max_count; ++k) {
        unsigned int idx  = feasible_indices(k);
        double       step = numerator(idx) * inv_diagH(idx);

        if(step < lb_u(idx)) step = lb_u(idx);
        if(step > ub_u(idx)) step = ub_u(idx);

        if(std::abs(step) > max_step_abs) {
          max_step_abs = std::abs(step);
          best_idx     = idx;
        }
      }
    }

    double lambda_opt_unclipped = numerator(best_idx) * inv_diagH(best_idx);
    double u_old                = u(best_idx);
    double u_new                = u_old + lambda_opt_unclipped;
    double delta;

    if (u_new < lb(best_idx)) {
      delta = lb(best_idx) - u_old;
      u(best_idx)    = lb(best_idx);
      lb_u(best_idx) = 0;
      ub_u(best_idx) = ub(best_idx) - lb(best_idx);
    } else if (u_new > ub(best_idx)) {
      delta = ub(best_idx) - u_old;
      u(best_idx)    = ub(best_idx);
      lb_u(best_idx) = lb(best_idx) - ub(best_idx);
      ub_u(best_idx) = 0;
    } else {
      delta           = lambda_opt_unclipped;
      u(best_idx)     = u_new;
      lb_u(best_idx) -= delta;
      ub_u(best_idx) -= delta;
    }
    numerator -= H.col(best_idx) * delta;
  }

  double obj_val = 0.5 * dot(u, q - numerator) - dot(q, u);

  return Rcpp::List::create(
    Rcpp::Named("x") = std::move(u),
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("objective.value") = obj_val
  );
}

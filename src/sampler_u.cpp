#include <Rcpp.h>
using namespace Rcpp;

namespace {

constexpr int kNumGrid = 250;

NumericVector build_continuous_grid(const double min_y, const double max_y, const double eps) {
  NumericVector grid(kNumGrid);
  const double start = min_y + eps;
  const double end = max_y - eps;
  const double step = (end - start) / (kNumGrid - 1.0);

  for (int j = 0; j < kNumGrid; ++j) {
    grid[j] = start + j * step;
  }

  return grid;
}

NumericMatrix build_basis_matrix(const NumericVector& theta,
                                 const NumericVector& h,
                                 const NumericVector& locations) {
  const int n = theta.size();
  const int m = locations.size();
  NumericMatrix basis(n, m);

  for (int i = 0; i < n; ++i) {
    const double h_i = h[i];
    const double theta_i = theta[i];
    for (int j = 0; j < m; ++j) {
      basis(i, j) = std::exp(theta_i * locations[j] + h_i);
    }
  }

  return basis;
}

NumericVector build_sum_vector(const NumericVector& u,
                               const NumericMatrix& basis) {
  const int n = u.size();
  const int m = basis.ncol();
  NumericVector sums(m);

  for (int j = 0; j < m; ++j) {
    double total = 0.0;
    for (int i = 0; i < n; ++i) {
      total += u[i] * basis(i, j);
    }
    sums[j] = total;
  }

  return sums;
}

double log_post_from_sums(const NumericVector& cont_sums,
                          const NumericVector& disc_sums,
                          const NumericVector& nstar,
                          const double alpha,
                          const double diff_z) {
  double integral_continuous = 0.0;
  for (int j = 0; j < cont_sums.size(); ++j) {
    integral_continuous += std::log1p(cont_sums[j]);
  }
  integral_continuous *= diff_z;

  double sum_discrete = 0.0;
  for (int j = 0; j < disc_sums.size(); ++j) {
    sum_discrete += nstar[j] * std::log1p(disc_sums[j]);
  }

  return -(alpha * integral_continuous + sum_discrete);
}

}  // namespace

// [[Rcpp::export]]
double log_post_u_cpp(const NumericVector& u,
                      const NumericVector& zstar,
                      const NumericVector& nstar,
                      const NumericVector& theta,
                      const double alpha,
                      const double min_y,
                      const double max_y,
                      const double eps,
                      const NumericVector& h) {
  const NumericVector cont_grid = build_continuous_grid(min_y, max_y, eps);;
  const double diff_z = (max_y - min_y - 2.0 * eps) / (kNumGrid - 1.0);
  const NumericMatrix cont_basis = build_basis_matrix(theta, h, cont_grid);
  const NumericMatrix disc_basis = build_basis_matrix(theta, h, zstar);
  const NumericVector cont_sums = build_sum_vector(u, cont_basis);
  const NumericVector disc_sums = build_sum_vector(u, disc_basis);

  return log_post_from_sums(cont_sums, disc_sums, nstar, alpha, diff_z);
}

// [[Rcpp::export]]
NumericVector sampler_u_cpp(NumericVector u,
                            const NumericVector& zstar,
                            const NumericVector& nstar,
                            const NumericVector& theta,
                            const double alpha,
                            const double delta,
                            const double min_y,
                            const double max_y,
                            const double eps,
                            const NumericVector& h) {
  const int n = u.size();
  const NumericVector cont_grid = build_continuous_grid(min_y, max_y, eps);
  const double diff_z = (max_y - min_y - 2.0 * eps) / (kNumGrid - 1.0);
  const NumericMatrix cont_basis = build_basis_matrix(theta, h, cont_grid);
  const NumericMatrix disc_basis = build_basis_matrix(theta, h, zstar);

  NumericVector cont_sums = build_sum_vector(u, cont_basis);
  NumericVector disc_sums = build_sum_vector(u, disc_basis);
  double current_log_post = log_post_from_sums(cont_sums, disc_sums, nstar, alpha, diff_z);

  for (int i = 0; i < n; ++i) {
    const double current_u_i = u[i];
    const double proposed_u_i = R::rgamma(delta, current_u_i / delta);
    const double delta_u = proposed_u_i - current_u_i;

    const double log_q_ratio =
      R::dgamma(current_u_i, delta, proposed_u_i / delta, true) -
      R::dgamma(proposed_u_i, delta, current_u_i / delta, true);

    NumericVector proposed_cont_sums = clone(cont_sums);
    for (int j = 0; j < proposed_cont_sums.size(); ++j) {
      proposed_cont_sums[j] += delta_u * cont_basis(i, j);
    }

    NumericVector proposed_disc_sums = clone(disc_sums);
    for (int j = 0; j < proposed_disc_sums.size(); ++j) {
      proposed_disc_sums[j] += delta_u * disc_basis(i, j);
    }

    const double proposed_log_post =
      log_post_from_sums(proposed_cont_sums, proposed_disc_sums, nstar, alpha, diff_z);

    const double log_ratio = proposed_log_post - current_log_post + log_q_ratio;

    if (std::log(R::runif(0.0, 1.0)) < log_ratio) {
      u[i] = proposed_u_i;
      cont_sums = proposed_cont_sums;
      disc_sums = proposed_disc_sums;
      current_log_post = proposed_log_post;
    }
  }

  return u;
}

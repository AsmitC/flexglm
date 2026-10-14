#include <Rcpp.h>
#include <vector>
#include <utility>
#include <limits>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

namespace {

constexpr int kNumGrid = 250;

NumericVector build_continuous_grid(const double min_y, const double max_y,
                                    const double eps) {
  NumericVector grid(kNumGrid);
  const double start = min_y + eps;
  const double end = max_y - eps;
  const double step = (end - start) / (kNumGrid - 1.0);

  for (int j = 0; j < kNumGrid; ++j) {
    grid[j] = start + j * step;
  }

  return grid;
}

// Store basis rows contiguously for componentwise proposals.
std::vector<double> build_basis_matrix(const NumericVector& theta,
                                       const NumericVector& h,
                                       const NumericVector& locations) {
  const int n = theta.size();
  const int m = locations.size();
  std::vector<double> basis(static_cast<size_t>(n) * m);

  for (int i = 0; i < n; ++i) {
    const double h_i = h[i];
    const double theta_i = theta[i];
    for (int j = 0; j < m; ++j) {
      basis[static_cast<size_t>(i) * m + j] =
        std::exp(theta_i * locations[j] + h_i);
    }
  }

  return basis;
}

NumericVector build_sum_vector(const NumericVector& u,
                               const std::vector<double>& basis,
                               const int m) {
  const int n = u.size();
  NumericVector sums(m);

  // Each sum still accumulates observations in the original order.
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      sums[j] += u[i] * basis[static_cast<size_t>(i) * m + j];
    }
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

// For a,b >= 0, integrating 1/(1+t) between a and b gives
// min(d/(1+a), d/(1+b)) <= log1p(b)-log1p(a)
//                       <= max(d/(1+a), d/(1+b)), where d=b-a.
// These bounds need no logarithms and become tight for coordinate updates in
// larger samples. Use the rounded proposed sum, just as the reference does.
inline void bound_log_difference(const double a, const double b,
                                 const double weight,
                                 double& lower, double& upper,
                                 double& max_sum, bool& safe) {
  const double d = b - a;
  const double at_a = d / (1.0 + a);
  const double at_b = d / (1.0 + b);
  lower += weight * std::min(at_a, at_b);
  upper += weight * std::max(at_a, at_b);
  max_sum = std::max(max_sum, std::max(a, b));
  safe = safe && a >= 0.0 && b >= 0.0;
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
  const NumericVector cont_grid = build_continuous_grid(min_y, max_y, eps);
  const double diff_z = (max_y - min_y - 2.0 * eps) / (kNumGrid - 1.0);
  const std::vector<double> cont_basis =
    build_basis_matrix(theta, h, cont_grid);
  const std::vector<double> disc_basis = build_basis_matrix(theta, h, zstar);
  const NumericVector cont_sums =
    build_sum_vector(u, cont_basis, cont_grid.size());
  const NumericVector disc_sums = build_sum_vector(u, disc_basis, zstar.size());

  return log_post_from_sums(cont_sums, disc_sums, nstar, alpha, diff_z);
}

// Updates u in place, including when a later CRM proposal is rejected.
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
  const std::vector<double> cont_basis =
    build_basis_matrix(theta, h, cont_grid);
  const std::vector<double> disc_basis = build_basis_matrix(theta, h, zstar);

  NumericVector cont_sums = build_sum_vector(u, cont_basis, cont_grid.size());
  NumericVector disc_sums = build_sum_vector(u, disc_basis, zstar.size());
  double current_log_post = 0.0;
  bool current_log_post_valid = false;
  double total_weight = alpha * diff_z * cont_sums.size();
  bool nonnegative_weights = alpha >= 0.0 && diff_z >= 0.0;
  for (int j = 0; j < nstar.size(); ++j) {
    total_weight += nstar[j];
    nonnegative_weights = nonnegative_weights && nstar[j] >= 0.0;
  }

  // Proposal sums
  NumericVector proposed_cont_sums(cont_sums.size());
  NumericVector proposed_disc_sums(disc_sums.size());
  for (int ii = 0; ii < 10; ++ii) {
    for (int i = 0; i < n; ++i) {
      const double current_u_i = u[i];
      const double proposed_u_i = R::rgamma(delta, current_u_i / delta);
      const double delta_u = proposed_u_i - current_u_i;

      const double log_q_ratio =
        R::dgamma(current_u_i, delta, proposed_u_i / delta, true) -
        R::dgamma(proposed_u_i, delta, current_u_i / delta, true);

      double cont_lower = 0.0, cont_upper = 0.0;
      double disc_lower = 0.0, disc_upper = 0.0;
      double max_sum = 0.0;
      bool safe = nonnegative_weights;
      for (int j = 0; j < proposed_cont_sums.size(); ++j) {
        proposed_cont_sums[j] = cont_sums[j] +
          delta_u * cont_basis[static_cast<size_t>(i) * cont_sums.size() + j];
        bound_log_difference(cont_sums[j], proposed_cont_sums[j], 1.0,
                             cont_lower, cont_upper, max_sum, safe);
      }

      for (int j = 0; j < proposed_disc_sums.size(); ++j) {
        proposed_disc_sums[j] = disc_sums[j] +
          delta_u * disc_basis[static_cast<size_t>(i) * disc_sums.size() + j];
        bound_log_difference(disc_sums[j], proposed_disc_sums[j], nstar[j],
                             disc_lower, disc_upper, max_sum, safe);
      }

      // MH step
      const double log_uniform = std::log(R::runif(0.0, 1.0));
      const double lower = alpha * (diff_z * cont_lower) + disc_lower;
      const double upper = alpha * (diff_z * cont_upper) + disc_upper;
      const double kernel_bound = total_weight * std::log1p(max_sum);
      // Allow for rounding error in the bounds
      const double margin = 64.0 * std::numeric_limits<double>::epsilon() *
        (cont_sums.size() + disc_sums.size() + 1.0) *
        (1.0 + 2.0 * kernel_bound + std::abs(log_q_ratio) +
         std::abs(log_uniform) + std::abs(lower) + std::abs(upper));
      safe = safe && std::isfinite(margin) && std::isfinite(lower) &&
        std::isfinite(upper);
      bool accept;
      bool proposed_log_post_valid = false;
      double proposed_log_post = 0.0;
      if (safe && log_uniform < log_q_ratio - upper - margin) {
        accept = true;
      } else if (safe && log_uniform >= log_q_ratio - lower + margin) {
        accept = false;
      } else {
        if (!current_log_post_valid) {
          current_log_post =
            log_post_from_sums(cont_sums, disc_sums, nstar, alpha, diff_z);
          current_log_post_valid = true;
        }
        proposed_log_post = log_post_from_sums(
          proposed_cont_sums, proposed_disc_sums, nstar, alpha, diff_z);
        proposed_log_post_valid = true;
        const double log_ratio =
          proposed_log_post - current_log_post + log_q_ratio;
        accept = log_uniform < log_ratio;
      }

      if (accept) {
        u[i] = proposed_u_i;
        std::swap(cont_sums, proposed_cont_sums);
        std::swap(disc_sums, proposed_disc_sums);
        current_log_post = proposed_log_post;
        current_log_post_valid = proposed_log_post_valid;
      }
    }
  }

  return u;
}

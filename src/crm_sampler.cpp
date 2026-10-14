#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

namespace {

// Jump-size grid
const std::vector<double>& get_s_grid(const int n_grid, const double eps) {
  static std::vector<double> grid;
  static int cached_n_grid = 0;
  static double cached_eps = 0.0;
  if (cached_n_grid != n_grid || cached_eps != eps) {
    grid.resize(n_grid);
    const double s_from = std::exp(-eps);
    const double s_to = std::exp(-5e-4);
    const double s_step = (s_to - s_from) / (n_grid - 1.0);
    for (int j = 0; j < n_grid; ++j) {
      grid[j] = -std::log(s_from + j * s_step);
    }
    cached_n_grid = n_grid;
    cached_eps = eps;
  }
  return grid;
}

inline double exp_sum_terms(const double* log_u_plus_h,
                            const double* theta,
                            const int n,
                            const double v) {
  double max_val = R_NegInf;

  for (int i = 0; i < n; ++i) {
    const double candidate = log_u_plus_h[i] + theta[i] * v;
    if (candidate > max_val) {
      max_val = candidate;
    }
  }

  double sum_exp = 0.0;
  for (int i = 0; i < n; ++i) {
    sum_exp += std::exp(log_u_plus_h[i] + theta[i] * v - max_val);
  }

  return std::exp(max_val) * sum_exp;
}

inline double compute_kernel_sum(const double s_j,
                                 const double* one_plus_psi_z,
                                 const int n_grid) {
  double kernel_sum = 0.0;
  for (int k = 0; k < n_grid; ++k) {
    kernel_sum += std::exp(-one_plus_psi_z[k] * s_j);
  }
  return kernel_sum;
}

NumericVector draw_RJ(const int M,
                      const std::vector<double>& s,
                      const std::vector<double>& Nv) {
  const int N = s.size();
  NumericVector xi(M);
  double xi_running = 0.0;
  for (int i = 0; i < M; ++i) {
    xi_running += R::rexp(1.0);
    xi[i] = xi_running;
  }

  NumericVector RJ(M);
  int iNv = N - 2;
  for (int i = 0; i < M; ++i) {
    while (iNv > 0 && Nv[iNv] < xi[i]) {
      --iNv;
    }
    RJ[i] = s[iNv + 1];
  }

  return RJ;
}

// Thinned Poisson process
NumericVector draw_RJ_thinned(const int M,
                              const std::vector<double>& s,
                              const std::vector<double>& bound_tail,
                              const std::vector<double>& endpoint_weights,
                              const std::vector<double>& rates,
                              const double min_rate) {
  // Remaining jumps are assigned to the grid floor
  NumericVector RJ(M, s[1]);
  int bin = s.size() - 2;
  int accepted = 0;
  double arrival = 0.0;
  while (accepted < M) {
    arrival += R::rexp(1.0);
    if (arrival > bound_tail[1]) break;
    while (bin > 1 && bound_tail[bin] < arrival) --bin;
    const double left = endpoint_weights[bin];
    const double right = endpoint_weights[bin + 1];
    const double endpoint = R::runif(0.0, 1.0) * (left + right) < left
      ? s[bin] : s[bin + 1];
    const int location = static_cast<int>(R::runif(0.0, 1.0) * rates.size());
    if (R::runif(0.0, 1.0) <
        std::exp(-(rates[location] - min_rate) * endpoint)) {
      RJ[accepted++] = s[bin + 1];
    }
  }
  return RJ;
}

NumericVector draw_RL(const NumericVector& RJ,
                      const std::vector<double>& z,
                      const std::vector<double>& one_plus_psi_z) {
  const int M = RJ.size();
  const int R = z.size();
  NumericVector RL(M);
  std::vector<double> cumulative_weights(R);
  double previous_rj = -1.0;
  double total_weight = 0.0;

  for (int m = 0; m < M; ++m) {
    const double rj_m = RJ[m];
    // Reuse location weights for repeated jumps
    if (rj_m != previous_rj) {
      total_weight = 0.0;
      for (int j = 0; j < R; ++j) {
        total_weight += std::exp(-one_plus_psi_z[j] * rj_m);
        cumulative_weights[j] = total_weight;
      }
      previous_rj = rj_m;
    }
    const double target = R::runif(0.0, 1.0) * total_weight;
    const auto selected = std::upper_bound(
      cumulative_weights.begin(), cumulative_weights.end(), target);
    RL[m] = z[selected == cumulative_weights.end() ? R - 1
                 : selected - cumulative_weights.begin()];
  }

  return RL;
}

NumericVector draw_Jstar(const NumericVector& nstar,
                         const NumericVector& psi_star) {
  const int n = nstar.size();
  NumericVector Jstar(n);
  for (int j = 0; j < n; ++j) {
    Jstar[j] = R::rgamma(nstar[j], 1.0 / (psi_star[j] + 1.0));
  }
  return Jstar;
}

List draw_crm_sample(const int M,
                     const std::vector<double>& s,
                     const std::vector<double>& Nv,
                     const std::vector<double>& z,
                     const std::vector<double>& one_plus_psi_z,
                     const NumericVector& zstar,
                     const NumericVector& nstar,
                     const NumericVector& psi_star,
                     const std::vector<double>& endpoint_weights,
                     const double min_rate,
                     const bool use_thinning) {
  const NumericVector RJ = use_thinning
    ? draw_RJ_thinned(M, s, Nv, endpoint_weights, one_plus_psi_z, min_rate)
    : draw_RJ(M, s, Nv);
  const NumericVector RL = draw_RL(RJ, z, one_plus_psi_z);
  const NumericVector Jstar = draw_Jstar(nstar, psi_star);

  return List::create(
    Named("RL") = RL,
    Named("RJ") = RJ,
    Named("zstar") = zstar,
    Named("Jstar") = Jstar
  );
}

}  // namespace

// [[Rcpp::export]]
List crm_sampler_cpp(const int M,
                     const NumericVector& u,
                     const NumericVector& zstar,
                     const NumericVector& nstar,
                     const NumericVector& theta,
                     const double alpha,
                     const double min_y,
                     const double max_y,
                     const double eps,
                     NumericVector h,
                     const int n_draws = 1,
                     const int n_grid = 3001,
                     const int r_grid = 3001) {
  if (M < 1) stop("M must be at least 1.");
  if (n_draws < 1) stop("n_draws must be at least 1.");
  if (n_grid < 2 || r_grid < 2) stop("Grid sizes must be at least 2.");
  if (!std::isfinite(min_y) || !std::isfinite(max_y) ||
      max_y <= min_y) {
    stop("The support must have finite, increasing bounds.");
  }
  if (!std::isfinite(alpha) || alpha <= 0.0) {
    stop("alpha must be positive and finite.");
  }
  if (u.size() == 0 || theta.size() != u.size()) {
    stop("theta must have the same positive length as u.");
  }
  if (zstar.size() != nstar.size()) {
    stop("zstar and nstar must have the same length.");
  }
  if (h.size() == 1 && u.size() > 1) h = rep(h[0], u.size());
  if (h.size() != u.size()) {
    stop("h must have length 1 or length equal to length(u).");
  }
  for (int i = 0; i < u.size(); ++i) {
    if (!std::isfinite(u[i]) || u[i] <= 0.0 ||
        !std::isfinite(theta[i]) || !std::isfinite(h[i])) {
      stop("u must be positive and finite; theta and h must be finite.");
    }
  }

  const std::vector<double>& s = get_s_grid(n_grid, eps);

  std::vector<double> z(r_grid);
  const double z_step = (max_y - min_y - 2.0 * eps) / (r_grid - 1.0);
  for (int j = 0; j < r_grid; ++j) {
    z[j] = min_y + eps + j * z_step;
  }

  const int n_theta = u.size();
  std::vector<double> log_u_plus_h(n_theta);
  for (int i = 0; i < n_theta; ++i) {
    log_u_plus_h[i] = std::log(u[i]) + h[i];
  }
  const double* log_u_plus_h_ptr = log_u_plus_h.data();
  const double* theta_ptr = REAL(theta);

  std::vector<double> one_plus_psi_z(r_grid);
  for (int j = 0; j < r_grid; ++j) {
    one_plus_psi_z[j] =
      1.0 + exp_sum_terms(log_u_plus_h_ptr, theta_ptr, n_theta, z[j]);
  }

  // Integration over the response support
  const double dz = z_step;
  std::vector<double> Nv(n_grid);
  const double min_rate =
    *std::min_element(one_plus_psi_z.begin(), one_plus_psi_z.end());
  std::vector<double> endpoint_weights(n_grid);
  for (int j = 0; j < n_grid; ++j) {
    endpoint_weights[j] = std::exp(-min_rate * s[j]) / s[j];
  }
  double tail_sum = 0.0;
  Nv[n_grid - 1] = 0.0;
  for (int j = n_grid - 2; j >= 0; --j) {
    tail_sum += (s[j + 1] - s[j]) *
      ((dz * r_grid * endpoint_weights[j] +
        dz * r_grid * endpoint_weights[j + 1]) / 2.0);
    Nv[j] = tail_sum * alpha;
  }
  // Use inversion when the expected number of candidates is large
  const bool use_thinning = std::isfinite(min_rate) &&
    std::isfinite(Nv[1]) && Nv[1] >= 0.0 && Nv[1] * n_draws <= 10000.0;
  if (!use_thinning) {
    tail_sum = 0.0;
    double fn_next = (dz / s[n_grid - 1]) *
      compute_kernel_sum(s[n_grid - 1], one_plus_psi_z.data(), r_grid);
    for (int j = n_grid - 2; j >= 0; --j) {
      const double fn_curr = (dz / s[j]) *
        compute_kernel_sum(s[j], one_plus_psi_z.data(), r_grid);
      tail_sum += (s[j + 1] - s[j]) * ((fn_curr + fn_next) / 2.0);
      Nv[j] = tail_sum * alpha;
      fn_next = fn_curr;
    }
  }

  NumericVector psi_star(zstar.size());
  for (int j = 0; j < zstar.size(); ++j) {
    psi_star[j] = exp_sum_terms(log_u_plus_h_ptr, theta_ptr, n_theta, zstar[j]);
  }

  if (n_draws == 1) {
    return draw_crm_sample(M, s, Nv, z, one_plus_psi_z, zstar, nstar, psi_star,
                           endpoint_weights, min_rate, use_thinning);
  }

  List draws(n_draws);
  for (int draw_idx = 0; draw_idx < n_draws; ++draw_idx) {
    draws[draw_idx] = draw_crm_sample(
      M, s, Nv, z, one_plus_psi_z, zstar, nstar, psi_star,
      endpoint_weights, min_rate, use_thinning
    );
  }
  return draws;
}

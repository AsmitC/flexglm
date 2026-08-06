#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstddef>

using namespace Rcpp;

namespace {

std::vector<double> make_s_grid(const double eps,
                                const int n_grid) {
  std::vector<double> s(n_grid);

  const double from = std::exp(-eps);
  const double to   = std::exp(-5e-4);
  const double step = (to - from) / (n_grid - 1.0);

  for (int j = 0; j < n_grid; ++j) {
    s[j] = -std::log(from + j * step);
  }

  return s;
}

std::vector<double> make_z_grid(const double min_y,
                                const double max_y,
                                const double eps,
                                const int r_grid) {
  std::vector<double> z(r_grid);

  const double step =
    (max_y - min_y - 2.0 * eps) / (r_grid - 1.0);

  for (int j = 0; j < r_grid; ++j) {
    z[j] = min_y + eps + j * step;
  }

  return z;
}

inline double log_sum_exp_terms(const std::vector<double>& log_u_plus_h,
                                const NumericVector& theta,
                                const double location) {
  const int n = theta.size();

  double max_val = R_NegInf;

  for (int i = 0; i < n; ++i) {
    const double val = log_u_plus_h[i] + theta[i] * location;
    if (val > max_val) max_val = val;
  }

  double sum_exp = 0.0;

  for (int i = 0; i < n; ++i) {
    sum_exp += std::exp(log_u_plus_h[i] + theta[i] * location - max_val);
  }

  return max_val + std::log(sum_exp);
}

inline double psi_at_location(const std::vector<double>& log_u_plus_h,
                              const NumericVector& theta,
                              const double location) {
  return std::exp(log_sum_exp_terms(log_u_plus_h, theta, location));
}

inline double lambda_at_location(
    const std::vector<double>& log_u_plus_h,
    const NumericVector& theta,
    const double location) {

  // lambda(z) = 1 + psi(z)
  return 1.0 + psi_at_location(log_u_plus_h, theta, location);
}


// Helper to compute sign of lambda'(z)
// Used to find min_z {lambda_z}
inline double scaled_lambda_derivative(
    const std::vector<double>& log_u_plus_h,
    const NumericVector& theta,
    const double location) {

  const int n = theta.size();
  double max_val = R_NegInf;

  for (int i = 0; i < n; ++i) {
    const double val = log_u_plus_h[i] + theta[i] * location;

    if (val > max_val) {
      max_val = val;
    }
  }

  double derivative_scaled = 0.0;

  for (int i = 0; i < n; ++i) {
    const double val = log_u_plus_h[i] + theta[i] * location;
    derivative_scaled += theta[i] * std::exp(val - max_val);
  }

  return derivative_scaled;
}

// Compute min_z {lambda_z}. Because lambda_z is convex,
// We only need to check the derivative + endpoints
double find_lambda_min(
    const std::vector<double>& log_u_plus_h,
    const NumericVector& theta,
    const double lower,
    const double upper) {

  const double derivative_lower =
    scaled_lambda_derivative(
      log_u_plus_h,
      theta,
      lower
    );

  const double derivative_upper =
    scaled_lambda_derivative(
      log_u_plus_h,
      theta,
      upper
    );

  double minimizer;

  if (derivative_lower >= 0.0) {
    // lambda_z increasing
    minimizer = lower;

  } else if (derivative_upper <= 0.0) {
    // lambda_z decreasing
    minimizer = upper;

  } else {
    // Interior root, use bisection
    double left = lower;
    double right = upper;

    for (int iteration = 0; iteration < 80; ++iteration) {
      const double midpoint = 0.5 * (left + right);
      const double derivative_midpoint =
        scaled_lambda_derivative(
          log_u_plus_h,
          theta,
          midpoint
        );

      if (derivative_midpoint > 0.0) {
        right = midpoint;
      } else {
        left = midpoint;
      }
    }
    minimizer = 0.5 * (left + right);
  }

  double lambda_min =
    lambda_at_location(
      log_u_plus_h,
      theta,
      minimizer
    );

  lambda_min = std::min(
    lambda_min,
    lambda_at_location(
      log_u_plus_h,
      theta,
      lower
    )
  );

  lambda_min = std::min(
    lambda_min,
    lambda_at_location(
      log_u_plus_h,
      theta,
      upper
    )
  );

  // Slightly decrease lambda_min found from optimization
  // to account for tolerance error
  return lambda_min * (1.0 - 1e-12);
}


// Compute C = E_G0[1 / {1 + psi(Z)}]
// where G0 ~ Unif[y_min - eps, y_max + eps]
double compute_C_uniform(
    const std::vector<double>& one_plus_psi_z,
    const double dz,
    const double support_width) {

  const int r_grid = static_cast<int>(one_plus_psi_z.size());

 // Approximate integral via trapezoidal rule
  double trapezoid_sum =
    0.5 / one_plus_psi_z.front() +
    0.5 / one_plus_psi_z.back();

  for (int j = 1; j < r_grid - 1; ++j) {
    trapezoid_sum += 1.0 / one_plus_psi_z[j];
  }

  const double integral = dz * trapezoid_sum;

  return integral / support_width;
}


// Generate S from the following model:
//
// Z ~ h_Z(z) \propto G_0(z) / (lambda_z)
// S | Z = z ~ Exp(rate = lambda_z)
// Accept S w.p. lambda_min / lambda_z
std::vector<double> sample_s_from_h(
    const int mc_size,
    const std::vector<double>& log_u_plus_h,
    const NumericVector& theta,
    const double lower,
    const double upper,
    const double lambda_min) {

  std::vector<double> s_mc;
  s_mc.reserve(mc_size);

  long long attempts = 0;

  while (
      static_cast<int>(s_mc.size()) < mc_size
  ) {

    // Sample z
    const double z_candidate = R::runif(lower, upper);

    // Compute lambda_z
    const double lambda_candidate =
      lambda_at_location(
        log_u_plus_h,
        theta,
        z_candidate
      );

    // Decide to accept / reject
    const double accept_probability = std::min(1.0, lambda_min / lambda_candidate);

    if (R::runif(0.0, 1.0) < accept_probability) {

      // R::rexp() takes scale
      const double s_draw = R::rexp(1.0 / lambda_candidate);

      if (std::isfinite(s_draw) && s_draw > 0.0) {
        s_mc.push_back(s_draw);
      }
    }

    ++attempts;

    if (attempts % 100000 == 0) {
      Rcpp::checkUserInterrupt();
    }
  }

  return s_mc;
}


// Approximate N(v) = alpha * C / mc_size * sum_m I(S_m > v) / S_m.
std::vector<double> estimate_Nv_mc(
    const std::vector<double>& v_grid,
    std::vector<double>& s_mc,
    const double alpha,
    const double C) {

  const int mc_size =
    static_cast<int>(s_mc.size());

  std::sort(
    s_mc.begin(),
    s_mc.end()
  );

  // suffix_sum[i] = sum_{m=i}^{mc_size-1} 1 / s_mc[m].
  std::vector<double> suffix_sum(
      mc_size + 1,
      0.0
  );

  for (int i = mc_size - 1; i >= 0; --i) {
    suffix_sum[i] = suffix_sum[i + 1] + 1.0 / s_mc[i];
  }

  std::vector<double> Nv(v_grid.size());

  for (std::size_t j = 0; j < v_grid.size(); ++j) {
    // Take first S satisfying S > v
    const auto first_greater =
      std::upper_bound(
        s_mc.begin(),
        s_mc.end(),
        v_grid[j]
      );

    const std::size_t index = static_cast<std::size_t>
      (first_greater - s_mc.begin());

    Nv[j] = alpha * C * suffix_sum[index] / static_cast<double>(mc_size);
  }

  return Nv;
}

NumericVector draw_RJ(const int M,
                      const std::vector<double>& s,
                      const std::vector<double>& Nv) {
  NumericVector RJ(M);

  const int n_grid = static_cast<int>(s.size());

  double xi = 0.0;
  int iNv = n_grid - 2;

  for (int i = 0; i < M; ++i) {
    xi += R::rexp(1.0);

    while (iNv > 0 && Nv[iNv] < xi) {
      --iNv;
    }

    RJ[i] = s[iNv + 1];
  }

  return RJ;
}

NumericVector draw_RL(
    const NumericVector& RJ,
    const std::vector<double>& z,
    const std::vector<double>& one_plus_psi_z) {

  const int M = RJ.size();
  const int r_grid = static_cast<int>(z.size());

  NumericVector RL(M);
  std::vector<double> weights(r_grid);

  for (int m = 0; m < M; ++m) {
    double total = 0.0;

    for (int j = 0; j < r_grid; ++j) {
      const double w =
        std::exp(-one_plus_psi_z[j] * RJ[m]);

      weights[j] = w;
      total += w;
    }

    const double target = R::runif(0.0, 1.0) * total;

    double cumulative = 0.0;
    int selected = r_grid - 1;

    for (int j = 0; j < r_grid; ++j) {
      cumulative += weights[j];

      if (cumulative > target) {
        selected = j;
        break;
      }
    }

    RL[m] = z[selected];
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

List draw_one_crm(
    const int M,
    const std::vector<double>& s,
    const std::vector<double>& Nv,
    const std::vector<double>& z,
    const std::vector<double>& one_plus_psi_z,
    const NumericVector& zstar,
    const NumericVector& nstar,
    const NumericVector& psi_star) {

  NumericVector RJ = draw_RJ(M, s, Nv);
  NumericVector RL = draw_RL(RJ, z, one_plus_psi_z);
  NumericVector Jstar = draw_Jstar(nstar, psi_star);

  return List::create(
    Named("RL") = RL,
    Named("RJ") = RJ,
    Named("zstar") = zstar,
    Named("Jstar") = Jstar
  );
}
} // namespace

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
                     const int n_grid = 1501,
                     const int r_grid = 1501,
                     const int mc_size = 1000) {
  const int n = u.size();

  if (M < 1) {
    stop("M must be at least 1.");
  }

  if (n_draws < 1) {
    stop("n_draws must be at least 1.");
  }

  if (n_grid < 2) {
    stop("n_grid must be at least 2.");
  }

  if (r_grid < 2) {
    stop("r_grid must be at least 2.");
  }

  if (mc_size < 1) {
    stop("mc_size must be at least 1.");
  }

  if (theta.size() != n) {
    stop("theta must have the same length as u.");
  }

  if (zstar.size() != nstar.size()) {
    stop("zstar and nstar must have the same length.");
  }

  if (h.size() == 1 && n > 1) {
    h = rep(h[0], n);
  }

  if (h.size() != n) {
    stop("h must have length 1 or length equal to length(u).");
  }

  if (min_y >= max_y) {
    stop("min_y must be less than max_y.");
  }

  if (eps <= 0) {
    stop("eps must be positive.");
  }

  for (int i = 0; i < n; ++i) {
    if (u[i] <= 0) {
      stop("All values of u must be positive.");
    }
  }

  std::vector<double> s = make_s_grid(eps, n_grid);

  std::vector<double> z = make_z_grid(min_y, max_y, eps, r_grid);

  const double dz = (max_y - min_y - 2.0 * eps) / (r_grid - 1.0);

  std::vector<double> log_u_plus_h(n);

  for (int i = 0; i < n; ++i) {
    log_u_plus_h[i] = std::log(u[i]) + h[i];
  }

  std::vector<double> psi_z(r_grid);
  std::vector<double> one_plus_psi_z(r_grid);

  for (int j = 0; j < r_grid; ++j) {
    psi_z[j] = psi_at_location(log_u_plus_h, theta, z[j]);
    one_plus_psi_z[j] = 1.0 + psi_z[j];
  }

  const double z_lower = z.front();
  const double z_upper = z.back();
  const double support_width = z_upper - z_lower;

    // Compute C
    compute_C_uniform(
      one_plus_psi_z,
      dz,
      support_width
    );

  if (!std::isfinite(C) || C <= 0.0) {
    stop("Failed to calculate a valid CRM normalizing constant C.");
  }


  // Compute min_z lambda_z
  const double lambda_min =
    find_lambda_min(
      log_u_plus_h,
      theta,
      z_lower,
      z_upper
    );

  if (!std::isfinite(lambda_min) || lambda_min <= 0.0) {
    stop("Failed to calculate a valid minimum of 1 + psi(z).");
  }

  // MC for S | Z
  std::vector<double> s_mc =
    sample_s_from_h(
      mc_size,
      log_u_plus_h,
      theta,
      z_lower,
      z_upper,
      lambda_min
    );


  // Approximate Nv
  std::vector<double> Nv =
    estimate_Nv_mc(
      s,
      s_mc,
      alpha,
      C
    );

  NumericVector psi_star(zstar.size());

  for (int j = 0; j < zstar.size(); ++j) {
    psi_star[j] = psi_at_location(log_u_plus_h, theta, zstar[j]);
  }

  if (n_draws == 1) {
    return draw_one_crm(
      M, s, Nv, z, one_plus_psi_z, zstar, nstar, psi_star
    );
  }

  List draws(n_draws);

  for (int d = 0; d < n_draws; ++d) {
    draws[d] = draw_one_crm(
      M, s, Nv, z, one_plus_psi_z, zstar, nstar, psi_star
    );
  }

  return draws;
}

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

using namespace Rcpp;

namespace {

inline double g_transform(const double mu, const double m, const double M) {
  return std::log(mu - m) - std::log(M - mu);
}

inline double clamp_interior(const double x, const double lower, const double upper) {
  if (x < lower) return lower;
  if (x > upper) return upper;
  return x;
}

void update_columns(const NumericVector& spt,
                    const NumericVector& f0,
                    const NumericVector& theta,
                    const std::vector<int>& cols,
                    NumericMatrix& fTilt,
                    NumericVector& bPrime,
                    NumericVector& bPrime2,
                    NumericVector& bPrimeErr,
                    const NumericVector& mu,
                    const bool logsumexp) {
  const int sptN = spt.size();

  for (int cc = 0; cc < static_cast<int>(cols.size()); ++cc) {
    const int i = cols[cc];
    const double theta_i = theta[i];

    double denom = 0.0;

    if (logsumexp) {
      double max_log = R_NegInf;

      for (int j = 0; j < sptN; ++j) {
        const double log_term = std::log(f0[j]) + spt[j] * theta_i;
        if (log_term > max_log) max_log = log_term;
      }

      for (int j = 0; j < sptN; ++j) {
        const double val = std::exp(std::log(f0[j]) + spt[j] * theta_i - max_log);
        fTilt(j, i) = val;
        denom += val;
      }
    } else {
      for (int j = 0; j < sptN; ++j) {
        const double val = f0[j] * std::exp(spt[j] * theta_i);
        fTilt(j, i) = val;
        denom += val;
      }
    }

    for (int j = 0; j < sptN; ++j) {
      fTilt(j, i) /= denom;
    }

    double mean_val = 0.0;
    for (int j = 0; j < sptN; ++j) {
      mean_val += spt[j] * fTilt(j, i);
    }

    double var_val = 0.0;
    for (int j = 0; j < sptN; ++j) {
      const double diff = spt[j] - mean_val;
      var_val += diff * diff * fTilt(j, i);
    }

    bPrime[i] = mean_val;
    bPrime2[i] = var_val;
    bPrimeErr[i] = mean_val - mu[i];
  }
}

std::vector<int> which_not_converged(const LogicalVector& conv,
                                     const LogicalVector& maxedOut) {
  std::vector<int> out;
  out.reserve(conv.size());

  for (int i = 0; i < conv.size(); ++i) {
    if (!conv[i] && !maxedOut[i]) out.push_back(i);
  }

  return out;
}

} // namespace

// [[Rcpp::export]]
List getTheta_cpp(const NumericVector& spt,
                  const NumericVector& f0,
                  const NumericVector& mu,
                  Nullable<NumericMatrix> sampprobs = R_NilValue,
                  Nullable<IntegerVector> ySptIndex = R_NilValue,
                  Nullable<NumericVector> thetaStart = R_NilValue,
                  const double eps = 1e-10,
                  const int maxiter = 100,
                  const int maxhalf = 20,
                  const double maxtheta = 500.0,
                  const bool logit = true,
                  const bool logsumexp = false) {
  const int sptN = spt.size();
  const int n = mu.size();

  if (f0.size() != sptN) {
    stop("spt and f0 must be vectors of equal length.");
  }

  for (int j = 0; j < sptN; ++j) {
    if (f0[j] < 0) {
      stop("f0 values cannot be negative.");
    }
  }

  const double m = min(spt);
  const double M = max(spt);

  if (min(mu) < m || max(mu) > M) {
    stop("mu starting values must lie within the range of spt.");
  }

  const double boundary_eps =
    std::max(1e-12 * (M - m), std::numeric_limits<double>::epsilon());
  const double safe_lower = m + boundary_eps;
  const double safe_upper = M - boundary_eps;
  const double var_eps = 1e-12;

  NumericVector theta(n);

  if (thetaStart.isNotNull()) {
    NumericVector theta_start_vec(thetaStart);

    if (theta_start_vec.size() != n) {
      stop("thetaStart must be a vector with length equal length(mu)");
    }

    theta = clone(theta_start_vec);
  } else {
    std::fill(theta.begin(), theta.end(), 0.0);
  }

  NumericVector gMu(n);
  for (int i = 0; i < n; ++i) {
    const double mu_safe = clamp_interior(mu[i], safe_lower, safe_upper);
    gMu[i] = g_transform(mu_safe, m, M);
  }

  NumericMatrix fTilt(sptN, n);
  NumericVector bPrime(n);
  NumericVector bPrime2(n);
  NumericVector bPrimeErr(n);

  std::vector<int> all_cols(n);
  for (int i = 0; i < n; ++i) all_cols[i] = i;

  update_columns(
    spt, f0, theta, all_cols,
    fTilt, bPrime, bPrime2, bPrimeErr,
    mu, logsumexp
  );

  NumericVector thetaOld(n, NA_REAL);
  NumericVector bPrimeErrOld(n, NA_REAL);
  LogicalVector conv(n, false);
  LogicalVector maxedOut(n, false);

  for (int i = 0; i < n; ++i) {
    conv[i] =
      (std::abs(bPrimeErr[i]) < eps) ||
      (theta[i] == maxtheta && bPrimeErr[i] < 0) ||
      (theta[i] == -maxtheta && bPrimeErr[i] > 0);
  }

  std::vector<int> active = which_not_converged(conv, maxedOut);
  int iter = 0;

  while (!active.empty() && iter < maxiter) {
    ++iter;

    for (int kk = 0; kk < static_cast<int>(active.size()); ++kk) {
      const int i = active[kk];

      bPrimeErrOld[i] = bPrimeErr[i];
      thetaOld[i] = theta[i];

      double theta_new = theta[i];

      if (logit) {
        const double bPrimeSafe = clamp_interior(bPrime[i], safe_lower, safe_upper);
        const double bPrime2Safe = std::max(bPrime2[i], var_eps);

        const double tPrime =
          (M - m) / ((bPrimeSafe - m) * (M - bPrimeSafe)) * bPrime2Safe;

        if (R_finite(tPrime) && std::abs(tPrime) > var_eps) {
          theta_new =
            theta[i] -
            (g_transform(bPrimeSafe, m, M) - gMu[i]) / tPrime;
        } else {
          theta_new =
            theta[i] -
            bPrimeErr[i] / std::max(bPrime2[i], var_eps);
        }
      } else {
        theta_new =
          theta[i] -
          bPrimeErr[i] / std::max(bPrime2[i], var_eps);
      }

      if (theta_new > maxtheta) theta_new = maxtheta;
      if (theta_new < -maxtheta) theta_new = -maxtheta;

      theta[i] = theta_new;
    }

    std::vector<int> ss = active;
    int nhalf = 0;

    while (!ss.empty() && nhalf < maxhalf) {
      update_columns(
        spt, f0, theta, ss,
        fTilt, bPrime, bPrime2, bPrimeErr,
        mu, logsumexp
      );

      std::vector<int> bad;
      bad.reserve(ss.size());

      for (int kk = 0; kk < static_cast<int>(ss.size()); ++kk) {
        const int i = ss[kk];

        if (std::abs(bPrimeErr[i]) > std::abs(bPrimeErrOld[i])) {
          bad.push_back(i);
        }
      }

      if (!bad.empty()) {
        ++nhalf;

        for (int kk = 0; kk < static_cast<int>(bad.size()); ++kk) {
          const int i = bad[kk];
          theta[i] = (theta[i] + thetaOld[i]) / 2.0;
        }
      }

      ss = bad;
    }

    for (int kk = 0; kk < static_cast<int>(ss.size()); ++kk) {
      const int i = ss[kk];

      maxedOut[i] = true;
      theta[i] = thetaOld[i];
    }

    for (int kk = 0; kk < static_cast<int>(active.size()); ++kk) {
      const int i = active[kk];

      conv[i] =
        (std::abs(bPrimeErr[i]) < eps) ||
        (theta[i] == maxtheta && bPrimeErr[i] < 0) ||
        (theta[i] == -maxtheta && bPrimeErr[i] > 0);
    }

    active = which_not_converged(conv, maxedOut);
  }

  NumericMatrix fTiltSW(sptN, n);
  NumericVector bPrimeSW(n);
  NumericVector bPrime2SW(n);

  if (sampprobs.isNotNull()) {
    NumericMatrix samp(sampprobs);

    if (samp.nrow() != n || samp.ncol() != sptN) {
      stop("sampprobs must have length(mu) rows and length(spt) columns.");
    }

    for (int i = 0; i < n; ++i) {
      double denom = 0.0;

      for (int j = 0; j < sptN; ++j) {
        fTiltSW(j, i) = fTilt(j, i) * samp(i, j);
        denom += fTiltSW(j, i);
      }

      for (int j = 0; j < sptN; ++j) {
        fTiltSW(j, i) /= denom;
      }
    }

    for (int i = 0; i < n; ++i) {
      double mean_val = 0.0;

      for (int j = 0; j < sptN; ++j) {
        mean_val += spt[j] * fTiltSW(j, i);
      }

      double var_val = 0.0;

      for (int j = 0; j < sptN; ++j) {
        const double diff = spt[j] - mean_val;
        var_val += diff * diff * fTiltSW(j, i);
      }

      bPrimeSW[i] = mean_val;
      bPrime2SW[i] = var_val;
    }
  } else {
    fTiltSW = clone(fTilt);
    bPrimeSW = clone(bPrime);
    bPrime2SW = clone(bPrime2);
  }

  double llik = NA_REAL;

  if (ySptIndex.isNotNull()) {
    IntegerVector y_idx(ySptIndex);

    if (y_idx.size() != n) {
      stop("ySptIndex must have length equal to length(mu).");
    }

    llik = 0.0;

    for (int i = 0; i < n; ++i) {
      const int j = y_idx[i] - 1;

      if (j < 0 || j >= sptN) {
        stop("ySptIndex contains an index outside the support range.");
      }

      llik += std::log(fTiltSW(j, i));
    }
  }

  return List::create(
    Named("theta") = theta,
    Named("fTilt") = fTilt,
    Named("bPrime") = bPrime,
    Named("bPrime2") = bPrime2,
    Named("fTiltSW") = fTiltSW,
    Named("bPrimeSW") = bPrimeSW,
    Named("bPrime2SW") = bPrime2SW,
    Named("llik") = llik,
    Named("conv") = conv,
    Named("iter") = iter
  );
}

## Required packages -----
## gldrm, mvtnorm, extraDistr

#' Stable Log-Sum-Exp
#'
#' Computes \eqn{\log \sum_j \exp(x_j)} using a numerically stable
#' log-sum-exp calculation.
#'
#' @param x Numeric vector.
#'
#' @return A numeric scalar.
#'
#' @keywords internal
logSumExp <- function(x) {
  i <- which.max(x)
  m <- x[i]
  lse <- log1p(sum(exp(x[-i] - m))) + m
  lse
}

#' Solve for Exponential-Tilting Parameters
#'
#' Computes the observation-specific tilting parameters, normalizing constants,
#' and second derivatives for the exponential-tilting model.
#'
#' @param spt Numeric vector giving the support points of the response.
#' @param f0 Numeric vector giving the baseline mass or density evaluated on
#' \code{spt}.
#' @param mu Numeric vector of target conditional means.
#' @param thtst Optional numeric vector of starting values for the tilting
#' parameters.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{bpr2}}{Second derivatives of the cumulant function.}
#'   \item{\code{tht}}{Tilting parameters.}
#'   \item{\code{btht}}{Normalizing constants \eqn{b(\theta_i)}.}
#' }
#'
#' @keywords internal
tht_sol <- function(spt, f0, mu, thtst) {
  out <- getTheta(
    spt = spt,
    f0 = f0,
    mu = mu,
    sampprobs = NULL,
    ySptIndex = NULL,
    thetaStart = thtst,
    thetaControl = gldrm::theta.control(logsumexp = TRUE)
  )
  tht <- out$theta
  bpr2 <- out$bPrime2
  lw <- sweep(tcrossprod(out$theta, spt), 2, log(f0), "+")
  btht <- apply(lw, 1, logSumExp)
  return(list(bpr2 = bpr2, tht = tht, btht = btht))
}

#' Compute Dirichlet Proposal Parameters
#'
#' Constructs the Dirichlet proposal parameter vector used in the update of the
#' baseline distribution.
#'
#' @param y Numeric response vector.
#' @param tht Numeric vector of current tilting parameters.
#' @param btht Numeric vector of current normalizing constants.
#' @param dir_pr_parm Numeric vector of Dirichlet prior parameters.
#' @param ind_mt Indicator matrix whose \eqn{(i,j)} entry identifies whether
#' response observation \eqn{i} equals support point \eqn{j}.
#'
#' @return Numeric vector of Dirichlet proposal parameters.
#'
#' @keywords internal
dir_parm <- function(y, tht, btht, dir_pr_parm, ind_mt) {
  wgt <- (1 / exp(tht * y - btht)) / sum(1 / exp(tht * y - btht))
  wgt <- as.numeric(wgt)
  parm <- dir_pr_parm + colSums(ind_mt * wgt)

  return(parm)
}

#' Function for log Dir pdf for x ~ Dir(a)
#'
#' @param x (vector) value at which to evaluate pdf
#' @param a (vector) Dirichlet parameter
#'
#' @return (vector) log Dirichlet pdf
#' @keywords internal
lDir <- function(x, a) {
  a0 <- sum(a)
  logB <- sum(lgamma(a)) - lgamma(a0)
  f <- sum((a - 1) * log(x)) - logB
  return(f)
}

#' Function for calculating \{f0(y_i), i = 1(1)n\}
#'
#' @param y (vector) response variable
#' @param spt (vector) support, \{s_j, j = 1(1)l\}, of response variable y
#' @param f0 (vector) \{f0(s_j), j=1(1)l\}: free parameter, f0(.) centering
#' density
#'
#' @return (vector) \{f0(y_i), i=1(1)n\}
#' @keywords internal
f0y <- function(y, spt, f0) {
  ind_mt <- outer(spt, y, `==`)
  f0_y <- colSums(f0 * ind_mt)
  return(f0_y)
}

#' Function for finding variance-covariance matrix of beta
#'
#' @param X design (matrix)
#' @param mu (vector) \{E(y_i | x_i), i=1(1)n\}
#' @param bpr2 (vector) \{second derivative of b(theta_i), i=1(1)n\}
#' @param rho MCMC update step size, a (scalar) in (0, 1]
#'
#' @return variance-covariance (matrix) of beta
#' @keywords internal
Sigma_beta <- function(X, mu, bpr2, rho, linkfun, mu.eta) {
  eta <- linkfun(mu)
  gpr <- as.numeric(1 / mu.eta(eta)) # Updated for general link
  gprsq_bpr2 <- gpr^2 * bpr2
  Xstar <- X / sqrt(gprsq_bpr2)
  info_mt <- t(Xstar) %*% Xstar
  Sigma <- rho^2 * solve(info_mt)

  return(Sigma)
}

#' Function for updating f0
#'
#' @param y response variable
#' @param spt y support
#' @param cr_f0 current f0 = \{f0(s_j), j=1(1)l\}
#' @param cr_f0y current f0(y) = \{f0(y_i), i=1(1)n\}
#' @param cr_dir_parm current dirichlet proposal parameter
#' @param cr_mu current mu
#' @param cr_tht current theta
#' @param cr_bpr2 current \{second derivative of b(theta_i), i=1(1)n\}
#' @param cr_btht current b(theta)
#' @param dir_pr_parm dirichlet prior parameter
#' @param ind_mt indicator matrix with (i,j)th element: 1(y_i == s_j)
#'
#' @importFrom extraDistr ddirichlet
#' @keywords internal
f0_update <- function(y,
                      spt,
                      cr_f0,
                      cr_f0y,
                      cr_dir_parm,
                      cr_mu,
                      cr_tht,
                      cr_bpr2,
                      cr_btht,
                      dir_pr_parm,
                      ind_mt) {
  l <- length(spt)
  eps <- .Machine$double.eps

  for (j in 1:l) {
    pr_f0 <- cr_f0
    pr_f0[j] <- rbeta(1, cr_dir_parm[j], sum(cr_dir_parm[-j]))
    pr_f0[-j] <- pr_f0[-j] * (1 - pr_f0[j]) / sum(pr_f0[-j])
    c1 <- cr_mu > spt[1] + eps
    c2 <- cr_mu < spt[l] - eps

    if (all(c1 & c2)) {
      out <- tht_sol(spt, pr_f0, cr_mu, cr_tht)
      pr_tht <- out$tht
      pr_btht <- out$btht
      pr_bpr2 <- out$bpr2
      pr_dir_parm <- dir_parm(y, pr_tht, pr_btht, dir_pr_parm, ind_mt)
      pr_f0y <- f0y(y, spt, pr_f0)

      pr_llik <- sum(pr_tht * y - pr_btht + log(pr_f0y))
      cr_llik <- sum(cr_tht * y - cr_btht + log(cr_f0y))

      R_pf0 <- sum((dir_pr_parm - 1) * log(pr_f0 / cr_f0))

      cr_qf0 <- lDir(cr_f0, pr_dir_parm)
      pr_qf0 <- lDir(pr_f0, cr_dir_parm)
      alp <- min(0, (pr_llik - cr_llik + R_pf0 + cr_qf0 - pr_qf0))

      if (log(runif(1)) < alp) {
        cr_f0 <- pr_f0
        cr_tht <- pr_tht
        cr_btht <- pr_btht
        cr_f0y <- pr_f0y
        cr_bpr2 <- pr_bpr2
        acc_f0 <- TRUE
      } else {
        acc_f0 <- FALSE
      }
    } else {
      acc_f0 <- FALSE
    } # Closes line 141
  } # Closes line 129
  return(
    list(
      cr_f0   = cr_f0,
      cr_f0y  = cr_f0y,
      cr_tht  = cr_tht,
      cr_btht = cr_btht,
      cr_bpr2 = cr_bpr2,
      acc_f0  = acc_f0
    )
  )
}


#' Function for updating beta (jointly)
#'
#' @param X design matrix
#' @param y response variable
#' @param spt support of y
#' @param cr_bt current beta
#' @param cr_Sig variance-covariance matrix for current beta
#' @param cr_f0 current f0 = \{f0(s_j), j=1(1)l\}
#' @param cr_tht current theta
#' @param cr_bpr2 current \{second derivative of b(theta_i), i=1(1)n\}
#' @param cr_btht current b(theta)
#' @param rho MCMC update step size, a (scalar) in (0, 1]
#'
#' @keywords internal
beta_update_joint <- function(X,
                              y,
                              spt,
                              cr_bt,
                              cr_Sig,
                              cr_f0,
                              cr_tht,
                              cr_bpr2,
                              cr_btht,
                              rho,
                              linkfun,
                              linkinv,
                              mu.eta,
                              mb,
                              Sb) {
  n <- dim(X)[1]
  l <- length(spt)

  pr_bt <- as.vector(mvtnorm::rmvnorm(1, mean = cr_bt, sigma = cr_Sig))
  pr_mu <- as.numeric(linkinv(X %*% pr_bt))

  if (
    all(is.finite(pr_mu)) &&
      sum(spt[1] <= pr_mu & pr_mu <= spt[l]) == n
  ) {
    out <- tht_sol(spt, cr_f0, pr_mu, cr_tht)
    pr_bpr2 <- out$bpr2
    pr_tht <- out$tht
    pr_btht <- out$btht
    pr_Sig <- Sigma_beta(X, pr_mu, pr_bpr2, rho, linkfun, mu.eta)

    pr_llik <- sum(pr_tht * y - pr_btht)
    cr_llik <- sum(cr_tht * y - cr_btht)
    pr_pbt <- mvtnorm::dmvnorm(pr_bt, mean = mb, sigma = Sb, log = TRUE)
    cr_pbt <- mvtnorm::dmvnorm(cr_bt, mean = mb, sigma = Sb, log = TRUE)
    cr_qbt <- mvtnorm::dmvnorm(cr_bt, mean = pr_bt, sigma = pr_Sig, log = TRUE)
    pr_qbt <- mvtnorm::dmvnorm(pr_bt, mean = cr_bt, sigma = cr_Sig, log = TRUE)

    alp <- min(0, pr_llik - cr_llik + pr_pbt - cr_pbt + cr_qbt - pr_qbt)

    if (log(runif(1)) < alp) {
      cr_bt <- pr_bt
      cr_tht <- pr_tht
      cr_btht <- pr_btht
      cr_bpr2 <- pr_bpr2
      cr_Sig <- pr_Sig
      acc_beta <- TRUE
    } else {
      acc_beta <- FALSE
    }
  } else {
    acc_beta <- FALSE
  }
  return(list(
    cr_bt    = cr_bt,
    cr_tht   = cr_tht,
    cr_btht  = cr_btht,
    cr_bpr2  = cr_bpr2,
    acc_beta = acc_beta
  ))
}


#' Function for updating beta (one at a time)
#'
#' @param X design matrix
#' @param y response variable
#' @param spt support of y
#' @param cr_bt current beta
#' @param cr_Sig variance-covariance matrix for current beta
#' @param cr_f0 current f0 = \{f0(s_j), j=1(1)l\}
#' @param cr_tht current theta
#' @param cr_bpr2 current \{second derivative of b(theta_i), i=1(1)n\}
#' @param cr_btht current b(theta)
#' @param rho MCMC update step size, a (scalar) in (0, 1]
#'
#' @keywords internal
beta_update_separate <- function(X,
                                 y,
                                 spt,
                                 cr_bt,
                                 cr_Sig,
                                 cr_f0,
                                 cr_tht,
                                 cr_bpr2,
                                 cr_btht,
                                 rho,
                                 linkfun,
                                 linkinv,
                                 mu.eta,
                                 mb,
                                 Sb) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  l <- length(spt)
  cr_sd <- sqrt(diag(cr_Sig))
  accepted <- logical(p)
  for (j in seq_len(p)) {
    pr_bt <- cr_bt
    pr_bt[j] <- cr_bt[j] + rnorm(1, mean = 0, sd = cr_sd[j])
    pr_mu <- as.numeric(linkinv(X %*% pr_bt)) # Updated for general link

    if (
      all(is.finite(pr_mu)) &&
        sum(spt[1] <= pr_mu & pr_mu <= spt[l]) == n
    ) {
      out <- tht_sol(spt, cr_f0, pr_mu, cr_tht)
      pr_tht <- out$tht
      pr_btht <- out$btht
      pr_bpr2 <- out$bpr2
      pr_Sig <- Sigma_beta(X, pr_mu, pr_bpr2, rho, linkfun, mu.eta)
      pr_sd <- sqrt(diag(pr_Sig))

      pr_llik <- sum(pr_tht * y - pr_btht)
      cr_llik <- sum(cr_tht * y - cr_btht)
      pr_pbt <- dnorm(pr_bt[j], mean = mb[j], sd = sqrt(Sb[j]), log = TRUE)
      cr_pbt <- dnorm(cr_bt[j], mean = mb[j], sd = sqrt(Sb[j]), log = TRUE)
      cr_qbt <- dnorm(cr_bt[j], mean = pr_bt[j], sd = pr_sd[j], log = TRUE)
      pr_qbt <- dnorm(pr_bt[j], mean = cr_bt[j], sd = cr_sd[j], log = TRUE)

      alp <- min(0, pr_llik - cr_llik + pr_pbt - cr_pbt + cr_qbt - pr_qbt)

      if (log(runif(1)) < alp) {
        cr_bt <- pr_bt
        cr_tht <- pr_tht
        cr_btht <- pr_btht
        cr_bpr2 <- pr_bpr2
        cr_Sig <- pr_Sig
        cr_sd <- pr_sd
        accepted[j] <- TRUE
      } else {
        accepted[j] <- FALSE
      }
    } else {
      accepted[j] <- FALSE
    }
  }
  return(list(
    cr_bt    = cr_bt,
    cr_tht   = cr_tht,
    cr_btht  = cr_btht,
    cr_bpr2  = cr_bpr2,
    acc_beta = mean(accepted)
  ))
}

### FUNCTIONS BELOW ARE FOR DPGLM ###
# -----------------------------------------------------------------------------#
#                               Update z                                      #
# -----------------------------------------------------------------------------#

z_sampler_unifK <- function(y, c0, crm.atoms, crm.jumps, tht, min_y, max_y,
                            eps) {
  n <- length(y)
  lower <- crm.atoms - c0
  upper <- crm.atoms + c0
  lower[lower < min_y] <- min_y + eps
  upper[upper > max_y] <- max_y - eps
  width <- upper - lower
  z <- numeric(n)
  for (i in 1:n) {
    indices <- y[i] >= lower & y[i] <= upper
    indx <- which(indices)
    log_prob <- -log(width[indx]) + (tht[i] * crm.atoms[indx]) +
      log(crm.jumps[indx])
    prob <- exp(log_prob - max(log_prob))
    if (sum(indices) == 1) {
      z[i] <- crm.atoms[indices]
    } else {
      z[i] <- sample(crm.atoms[indx], 1, prob = prob)
    }
  }
  return(z)
}

c0_silverman <- function(y) {
  n <- length(y)
  sd_y <- sd(y)
  iqr_y <- IQR(y)
  return(0.9 * min(sd_y, iqr_y / 1.34) * n^(-1 / 5))
}

#' Resample Unique Latent Atom Locations
#'
#' Collapses the current latent atom assignments into unique atom locations and
#' their multiplicities. These unique values are used as the fixed-location
#' component in the CRM update.
#'
#' @param z Numeric vector of latent atom assignments.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{zstar}}{Numeric vector of unique latent atom locations.}
#'   \item{\code{nstar}}{Numeric vector of multiplicities corresponding to
#'   \code{zstar}.}
#' }
#'
#' @keywords internal
resample_zstar <- function(z) {
  z_table <- table(z)
  zstar <- as.numeric(names(z_table))
  nstar <- as.numeric(z_table)
  return(list(zstar = zstar, nstar = nstar))
}

#' Log Posterior Kernel for Beta
#'
#' Evaluates the log posterior kernel for beta conditional on
#' the currently imputed values.
#'
#' @param beta Numeric vector of regression coefficients.
#' @param linkinv Inverse link function mapping linear predictors to means.
#' @param z Numeric vector of latent atom assignments.
#' @param X Numeric design matrix.
#' @param atoms Numeric vector of CRM atom locations.
#' @param jumps Numeric vector of CRM jump weights.
#' @param mu_beta Numeric vector giving the prior mean for \code{beta}.
#' @param sigma_beta Numeric vector giving the prior standard deviations for
#' \code{beta}.
#' @param h Numeric vector of log-likelihood contributions from dependence.
#' @param theta Optional tilting parameters at the supplied beta and CRM.
#'
#' @return Numeric scalar giving the log posterior kernel, up to
#' proportionality.
#'
#' @keywords internal
logpost_beta <- function(beta, linkinv, z, X, atoms, jumps, mu_beta,
                         sigma_beta, h, theta = NULL) {
  mu <- linkinv(X %*% beta)
  if (is.null(theta)) {
    theta <- getTheta(
      spt = atoms, f0 = jumps, mu = mu,
      ySptIndex = NULL, sampprobs = NULL
    )$theta
  }
  btheta <- b_theta(theta, atoms, jumps)
  # others const in beta
  log_post <- sum(theta * z + h - btheta) -
    sum((beta - mu_beta)^2 / (2 * sigma_beta^2))
  return(log_post)
}

#' Model Log Likelihood
#'
#' Computes the observed log likelihood contribution under
#' the currently imputed values.
#'
#' @param linkinv Inverse link function.
#' @param z Numeric vector of latent atom assignments.
#' @param X Numeric design matrix.
#' @param beta Numeric vector of regression coefficients.
#' @param atoms Numeric vector of CRM atom locations.
#' @param jumps Numeric vector of CRM jump weights.
#'
#' @return Numeric scalar giving the log likelihood.
#'
#' @keywords internal
loglik <- function(linkinv, z, X, beta, atoms, jumps) {
  n <- length(z)
  mu <- linkinv(X %*% beta)
  theta <- getTheta(
    spt = atoms, f0 = jumps, mu = mu, ySptIndex = NULL,
    sampprobs = NULL
  )$theta

  btheta <- b_theta(theta, atoms, jumps)

  f0_z <- numeric(n)
  for (i in 1:n) {
    f0_z[i] <- sum(jumps[z[i] == atoms])
  }
  loglik <- sum(theta * z - btheta + log(f0_z))
  return(loglik)
}

#' Log Posterior Kernel for u
#'
#' Evaluates the log posterior kernel for the auxiliary variables \code{u} used
#' in the CRM representation.
#'
#' @param u Numeric vector of auxiliary variables.
#' @param zstar Numeric vector of unique latent atom locations.
#' @param nstar Numeric vector of multiplicities corresponding to \code{zstar}.
#' @param theta Numeric vector of observation-specific tilting parameters.
#' @param alpha Numeric concentration parameter for the CRM prior.
#' @param min_y Numeric lower bound of the response support.
#' @param max_y Numeric upper bound of the response support.
#' @param eps Numeric boundary padding constant.
#' @param h Numeric vector of log-likelihood contributions from dependence.
#'
#' @return Numeric scalar giving the log posterior kernel for \code{u}
#' up to proportionality.
#'
#' @keywords internal
log_post_u <- function(u, zstar, nstar, theta, alpha, min_y, max_y, eps, h) {
  if (length(h) == 1L) h <- rep(h, length(u))
  log_post_u_cpp(u, zstar, nstar, theta, alpha, min_y, max_y, eps, h)
}

#' Update u
#'
#' Performs ten MH sweeps for the auxiliary variables \code{u}.
#' Each component is proposed from a Gamma random walk. The compiled sampler
#' reuses cached sums and bounds acceptance decisions before evaluating the
#' full log posterior when needed.
#'
#' @param u Numeric vector of current auxiliary variables.
#' @param zstar Numeric vector of unique latent atom locations.
#' @param nstar Numeric vector of multiplicities corresponding to \code{zstar}.
#' @param theta Numeric vector of observation-specific tilting parameters.
#' @param alpha Numeric concentration parameter for the CRM prior.
#' @param delta Numeric tuning parameter for the Gamma proposal.
#' @param min_y Numeric lower bound of the response support.
#' @param max_y Numeric upper bound of the response support.
#' @param eps Numeric boundary padding constant.
#' @param h Numeric vector of log-likelihood contributions from dependence.
#'
#' @return Numeric vector of updated auxiliary variables.
#'
#' @keywords internal
sampler_u <- function(u, zstar, nstar, theta, alpha, delta,
                      min_y, max_y, eps, h = 0) {
  if (length(h) == 1L) {
    h <- rep(h, length(u))
  }

  sampler_u_cpp(
    u = u,
    zstar = zstar,
    nstar = nstar,
    theta = theta,
    alpha = alpha,
    delta = delta,
    min_y = min_y,
    max_y = max_y,
    eps = eps,
    h = h
  )
}

#' Compute b(theta) for discrete base measure
#'
#' Calculates \eqn{b(\theta) = \log \left( \int \exp(\theta z) dG_0(z) \right)}
#' for a discrete base measure \eqn{G_0}
#' with support points \code{spt} and weights \code{f0}, using log-sum-exp trick
#' for numerical stability.
#'
#' @param theta Numeric vector of parameter values.
#' @param spt Numeric vector of support points for the base measure.
#' @param f0 Numeric vector of weights or density values at each support point.
#'
#' @return Numeric vector of \eqn{b(\theta)} values, one for each \code{theta}.
#' @export
b_theta <- function(theta, spt, f0) {
  log_f0_dz <- log(f0) # log density values

  # Create the log weights matrix of dimensions length(theta) x length(spt)
  # (i, j)th element is theta[i]*spt[j] + log(f0[j])
  log_weights <- outer(theta, spt, "*") +
    matrix(log_f0_dz, nrow = length(theta), ncol = length(spt), byrow = TRUE)

  # Use log-sum-exp trick to compute b(theta) for each theta
  result <- apply(log_weights, 1, function(x) {
    m <- max(x)
    m + log(sum(exp(x - m)))
  })

  return(as.vector(result)) # Ensure the output is a vector
}

#' Sample a CRM Proposal
#'
#' Generates a proposed CRM update using Poisson thinning, with grid integration
#' and tail inversion when the expected thinning work is large.
#'
#' @param M Integer. Number of random CRM atoms to generate.
#' @param u Numeric vector of auxiliary variables.
#' @param zstar Numeric vector of unique occupied latent atom locations.
#' @param nstar Numeric vector of multiplicities corresponding to \code{zstar}.
#' @param theta Numeric vector of observation-specific tilting parameters.
#' @param alpha Numeric concentration parameter for the CRM prior.
#' @param min_y Numeric lower bound of the response support.
#' @param max_y Numeric upper bound of the response support.
#' @param eps Numeric boundary padding constant.
#' @param h Numeric vector of log-likelihood contributions from dependence.
#'
#' @param n_draws Number of CRM draws sharing the same setup. Defaults to 1.
#' @param n_grid Number of points in the jump-size grid. Defaults to 3001.
#' @param r_grid Number of points in the response-support grid. Defaults to
#' 3001.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{RL}}{Random atom locations.}
#'   \item{\code{RJ}}{Random jump sizes corresponding to \code{RL}.}
#'   \item{\code{zstar}}{Fixed occupied atom locations.}
#'   \item{\code{Jstar}}{Jump sizes for fixed occupied atom locations.}
#' }
#'
#' @keywords internal
crm_sampler <- function(
  M,
  u,
  zstar,
  nstar,
  theta,
  alpha,
  min_y,
  max_y,
  eps,
  h = 0,
  n_draws = 1L,
  n_grid = 3001L,
  r_grid = 3001L
) {
  crm_sampler_cpp(
    M = M,
    u = u,
    zstar = zstar,
    nstar = nstar,
    theta = theta,
    alpha = alpha,
    min_y = min_y,
    max_y = max_y,
    eps = eps,
    h = h,
    n_draws = n_draws,
    n_grid = n_grid,
    r_grid = r_grid
  )
}

### FUNCTIONS BELOW ARE FOR CDPGLM ###

#' Compute observation-specific marginal CDFs
#'
#' Computes observation-specific marginal CDF values
#' \deqn{
#' F_{x_j}(y_j)
#' }
#' using a CRM-induced mixture of truncated uniform kernels.
#'
#' The mixture weights are covariate-dependent and vary across observations
#' through \eqn{\theta_j}.
#'
#' @param y Numeric vector of outcomes with length \eqn{n}.
#' @param crm.atoms Vector of CRM atoms \eqn{(z_{\ell})}.
#' @param crm.jumps Vector of CRM jump sizes \eqn{(J_{\ell})}.
#' @param theta Numeric vector of covariate effects with length \eqn{n}.
#' @param c0 Half-width of the uniform kernel.
#' @param min_y Global lower bound of the outcome support.
#' @param max_y Global upper bound of the outcome support.
#'
#' @return
#' Numeric vector of marginal CDF values
#' \deqn{
#' F_{x_j}(y_j).
#' }
#'
#' @keywords internal
marginal_cdf_unif_kde <- function(y, crm.atoms, crm.jumps, theta, c0,
                                  min_y, max_y) {
  eps <- 1e-6
  n <- length(y)

  ## 1. Construct truncated kernel bounds (shared across observations)
  lower <- crm.atoms - c0
  upper <- crm.atoms + c0

  lower[lower < min_y] <- min_y + eps
  upper[upper > max_y] <- max_y - eps

  width <- upper - lower

  ## 2. Evaluate CDF observation by observation
  Fy <- numeric(n)

  for (j in seq_len(n)) {
    ## Observation-specific mixture weights (log-scale)
    logw <- log(crm.jumps) + theta[j] * crm.atoms
    logw <- logw - max(logw) # log-sum-exp stabilization
    w <- exp(logw)
    w <- w / sum(w)

    ## Mixture CDF at y_j
    yj <- y[j]
    Fy[j] <- sum(
      w * ifelse(yj < lower, 0, ifelse(yj > upper, 1, (yj - lower) / width))
    )
  }

  Fy
}


#' Compute the log-density of a Gaussian copula
#'
#' Computes the log-density of a Gaussian copula with either
#' exchangeable or AR1 correlation structure.
#'
#' The copula density is evaluated at a vector of marginal probabilities
#' \deqn{
#' u = (u_1, \ldots, u_n),
#' }
#' where
#' \deqn{
#' u_j = F_x(y_{ij}).
#' }
#'
#' @param u Numeric vector of marginal CDF values in \eqn{(0,1)}.
#' @param rho Dependence parameter, interpreted as the pairwise correlation.
#' @param corr Type of correlation structure. Must be one of \code{ex} or
#' \code{ar1}.
#'
#' @return
#' Scalar value of
#' \deqn{
#' \log c_{\rho}(u).
#' }
#'
#' @keywords internal
log_copula_gaussian <- function(u, rho, corr) {
  n <- length(u)

  # avoid boundary
  eps <- 1e-10
  u <- pmin(pmax(u, eps), 1 - eps)

  ## Transform to latent Gaussian scale
  z <- qnorm(u)

  if (corr == "ex") {
    Sigma <- matrix(rho, n, n)
    diag(Sigma) <- 1
  } else {
    Sigma <- outer(1:n, 1:n, function(i, j) rho^abs(i - j))
  }

  Sigma_inv <- solve(Sigma)
  logdet <- determinant(Sigma, logarithm = TRUE)$modulus

  ## Gaussian copula log-density
  # log c(u) = -1/2 log|Sigma| - 1/2 z' (Sigma^{-1} - I) z
  val <- -0.5 * logdet -
    0.5 * t(z) %*% (Sigma_inv - diag(n)) %*% z

  as.numeric(val)
}

#' Compute observation-level Gaussian copula log-likelihood contributions
#'
#' Marginal CDFs
#' \deqn{
#' F_j = F_{x_j}(y_j)
#' }
#' are computed internally via a CRM-induced mixture of truncated
#' uniform kernels using \code{marginal_cdf_unif_kde()}.
#'
#' Observations are coupled within groups using a Gaussian copula with
#' exchangeable correlation parameter \eqn{\rho}.
#'
#' For a group \eqn{i} with \eqn{n_i} observations, each observation
#' \eqn{j} in that group receives an equal share of the group copula
#' log-density:
#'
#' \deqn{
#' \log c_{\rho}(F_{i1}, \ldots, F_{in_i}) / n_i
#' }
#'
#' @param y Numeric vector of outcomes.
#' @param group_index Vector indicating group membership for each   observation
#' in \code{y}.
#' @param rho Copula dependence parameter.
#' @param crm.atoms Vector of CRM atoms \eqn{(z_{\ell})}.
#' @param crm.jumps Vector of CRM jump sizes \eqn{(J_{\ell})}.
#' @param theta Numeric vector of covariate effects with length equal   to
#' \code{length(y)}.
#' @param c0 Half-width of the uniform kernel.
#' @param min_y Global lower bound of the outcome support.
#' @param max_y Global upper bound of the outcome support.
#'
#' @return
#' Numeric vector \code{logcop_obs} of observation-level copula
#' log-likelihood contributions with the same length as \code{y}.
#'
#' @keywords internal
log_copula_contribution_by_obs <- function(y, group_index, rho, corr,
                                           crm.atoms, crm.jumps, theta, c0,
                                           min_y, max_y) {
  ## 1. Compute marginal CDFs internally
  Fy <- marginal_cdf_unif_kde(
    y = y,
    crm.atoms = crm.atoms, crm.jumps = crm.jumps, theta = theta, c0 = c0,
    min_y = min_y, max_y = max_y
  )

  ## 2. Initialize observation-level output
  logcop_obs <- numeric(length(y))

  ## 3. Compute copula contribution per group, then distribute to observations
  groups <- unique(group_index)

  for (g in groups) {
    idx <- which(group_index == g)
    u <- Fy[idx]
    n_i <- length(u)

    logc_group <- log_copula_gaussian(u, rho, corr)

    # Assign equal share to each observation in the group
    logcop_obs[idx] <- logc_group / n_i
  }

  logcop_obs
}

# Bounds for the copula correlation
rho_bounds <- function(corr, group_index = NULL) {
  if (corr == "ar1") return(c(-1, 1))
  if (corr != "ex") stop("corr must be one of 'ex' or 'ar1'.")
  if (!length(group_index) || anyNA(group_index)) {
    stop("group_index must contain non-missing group labels.")
  }
  m <- max(table(group_index))
  c(if (m > 1L) -1 / (m - 1L) else -1, 1)
}

# Initial rho
check_rho_start <- function(rho, bounds, corr) {
  default <- mean(bounds)
  if (is.null(rho)) return(default)
  if (
    !is.numeric(rho) || length(rho) != 1L || !is.finite(rho) ||
      rho <= bounds[1] || rho >= bounds[2]
  ) {
    warning(sprintf(
      "rhoStart must lie in (%g, %g) for corr = '%s'; using %g.",
      bounds[1], bounds[2], corr, default
    ), call. = FALSE)
    return(default)
  }
  as.numeric(rho)
}

#' Log posterior kernel for rho
#'
#' Adds the copula contribution and a Beta prior on rho rescaled to (0, 1).
#'
#' @param rho Copula correlation, strictly inside \code{bounds}.
#' @param logc_obs Observation-specific log copula contributions.
#' @param shape1,shape2 Beta prior shape parameters.
#' @param bounds Lower and upper bounds for rho. Defaults to (0, 1).
#' @return Log posterior, up to an additive constant.
#' @keywords internal
log_post_rho <- function(rho, logc_obs, shape1, shape2, bounds = c(0, 1)) {
  if (!is.finite(rho) || rho <= bounds[1] || rho >= bounds[2]) return(-Inf)
  rho_unit <- (rho - bounds[1]) / diff(bounds)
  sum(logc_obs) + dbeta(rho_unit, shape1, shape2, log = TRUE) -
    log(diff(bounds))
}
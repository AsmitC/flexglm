#' Control Arguments for CDPGLM
#'
#' Creates a list of control parameters used by \code{\link{cdpglm}} and
#' \code{\link{cdpglmFit}}.
#'
#' @param burnin Number of burn-in MCMC iterations. Defaults to 100.
#' @param thin Factor by which to thin MCMC iterations. Defaults to 10.
#' @param save Number of MCMC samples to return. Defaults to 1000.
#' @param mb Prior mean for beta. Defaults to NULL.
#' @param Sb Prior scale / diagonal variance entries for beta. Defaults to NULL.
#' @param M Number of random CRM atoms proposed at each CRM update. Defaults to
#' 50.
#' @param alpha CRM mass parameter. Defaults to 1.
#' @param delta Tuning parameter for the \code{u} update. Defaults to 2.
#' @param c0 Bandwidth parameter used in the latent \code{z} update. Defaults to
#' NULL.
#' @param gamma Shrinkage parameter for the default prior variance on
#' \code{beta}. Defaults to 1. Ignored if \code{Sb} is specified.
#' @param mu0 Mean of the reference distribution. Defaults to \code{NULL}, in
#' which case \code{mean(y)} is used internally.
#' @param spt Theoretical support for the response in the DP-GLM. This should be
#' a length-2 vector giving the lower and upper bounds. Defaults to \code{NULL},
#' in which case the support is inferred internally from the data as
#' \code{c(min(y)-eps, max(y)+eps)}.
#' @param eps Padding used when constructing default support bounds. Defaults to
#' 1e-6.
#' @param n_grid Number of grid points for CRM jump sizes. Defaults to 3001.
#' @param r_grid Number of grid points for CRM locations. Defaults to 3001.
#' @param rhoStart Numeric. Starting value for the copula dependence parameter
#' \eqn{\rho}. Defaults to the midpoint of its admissible interval. Invalid
#' values produce a warning and are replaced by this default.
#' @param corr Character string specifying copula correlation structure. Options
#' are \code{"ex"} for exchangeable correlation and \code{"ar1"} for AR(1)
#' correlation.
#' @param rho_proposal_sd Numeric. Proposal standard deviation for the
#' random-walk update on the logit scale of the rescaled rho. Defaults to
#' \code{0.1}.
#' @param rho_prior_shape Numeric vector of length two. Parameters for the Beta
#' prior on rho rescaled to (0, 1). Defaults to \code{c(8, 2)}. The rho interval
#' is (-1, 1) for AR(1), and (-1/(m - 1), 1) for exchangeable correlation, where
#' m is the largest group size. If all groups are singletons, (-1, 1) is used.
#' @param betaStart Initial value for \code{beta}. Defaults to NULL.
#' @param varbetaStart Initial value for the proposal covariance of \code{beta}.
#' Defaults to NULL.
#' @param thetaStart Initial value for \code{theta}. Defaults to NULL.
#' @param crmStart Initial value for the CRM state, given as a list with
#' components \code{z.tld} and \code{J.tld}. Defaults to NULL.
#' @param seed Random seed. Defaults to NULL.
#' @param robust Provides numerical stability. Defaults to \code{TRUE}. See
#' details for more.
#'
#' @return An object of class \code{"cdpglmControl"}; a list of control
#' parameters passed to the internal CDPGLM fitting routine.
#'
#' @examples
#' ctrl <- cdpglm.control(
#'   burnin = 100,
#'   thin = 5,
#'   save = 500,
#'   spt = c(0, 1),
#'   corr = "ex",
#'   seed = 123
#' )
#'
#' @export
cdpglm.control <- function(burnin = 100, thin = 10, save = 1000,
                           mb = NULL, Sb = NULL,
                           M = 50, alpha = 1, delta = 2, c0 = NULL,
                           gamma = 1, mu0 = NULL, spt = NULL, eps = 1e-6,
                           n_grid = 3001L, r_grid = 3001L,
                           rhoStart = NULL, corr = "ex",
                           rho_proposal_sd = 0.1,
                           rho_prior_shape = c(8, 2),
                           betaStart = NULL, varbetaStart = NULL,
                           thetaStart = NULL, crmStart = NULL,
                           seed = NULL, robust = TRUE) {
  if (burnin < 0 || floor(burnin) != burnin) {
    stop("Number of burn-in samples must be an integer >= 0.", call. = FALSE)
  }
  if (thin < 1 || floor(thin) != thin) {
    stop("Thin must be an integer >= 1.", call. = FALSE)
  }
  if (save < 1 || floor(save) != save) {
    stop("Number of saved iterations must be an integer >= 1.", call. = FALSE)
  }

  if (!(corr %in% c("ex", "ar1"))) stop("corr must be one of c('ex', 'ar1')")

  if (
    !is.null(rhoStart) &&
      (!is.numeric(rhoStart) || length(rhoStart) != 1L ||
         !is.finite(rhoStart) || rhoStart <= -1 || rhoStart >= 1)
  ) {
    warning("Invalid rhoStart; using the correlation-specific default.",
            call. = FALSE)
    rhoStart <- NULL
  }

  if (
    !is.numeric(rho_proposal_sd) || length(rho_proposal_sd) != 1L ||
      !is.finite(rho_proposal_sd) || rho_proposal_sd <= 0
  ) {
    stop("rho_proposal_sd must be a positive numeric value.", call. = FALSE)
  }

  if (
    !is.numeric(rho_prior_shape) || length(rho_prior_shape) != 2L ||
      any(!is.finite(rho_prior_shape)) || any(rho_prior_shape <= 0)
  ) {
    stop(
      paste0("rho_prior_shape must be a numeric vector of length 2 ",
             "with positive entries."),
      call. = FALSE
    )
  }

  if (
    !is.numeric(n_grid) ||
      length(n_grid) != 1L ||
      !is.finite(n_grid) ||
      n_grid < 2 ||
      n_grid != floor(n_grid)
  ) {
    stop("n_grid must be an integer >= 2.")
  }

  if (
    !is.numeric(r_grid) ||
      length(r_grid) != 1L ||
      !is.finite(r_grid) ||
      r_grid < 2 ||
      r_grid != floor(r_grid)
  ) {
    stop("r_grid must be an integer >= 2.")
  }

  n_grid <- as.integer(n_grid)
  r_grid <- as.integer(r_grid)

  ctrl <- list(
    burnin = burnin,
    thin = thin,
    save = save,
    mb = mb,
    Sb = Sb,
    M = M,
    alpha = alpha,
    delta = delta,
    c0 = c0,
    gamma = gamma,
    mu0 = mu0,
    spt = spt,
    eps = eps,
    n_grid = n_grid,
    r_grid = r_grid,
    rhoStart = rhoStart,
    corr = corr,
    rho_proposal_sd = rho_proposal_sd,
    rho_prior_shape = rho_prior_shape,
    betaStart = betaStart,
    varbetaStart = varbetaStart,
    thetaStart = thetaStart,
    crmStart = crmStart,
    seed = seed,
    robust = robust
  )

  class(ctrl) <- "cdpglmControl"
  ctrl
}

#' Main MCMC Fitting Function for CDPGLM
#'
#' Performs MCMC to fit the CDPGLM model.
#'
#' @param formula Model formula supplied to \code{\link{cdpglm}}.
#' @param data Data frame containing the variables appearing in \code{formula}.
#' @param X Numeric design matrix.
#' @param y Numeric response vector.
#' @param group_index Vector identifying the dependence group for each
#' observation.
#' @param link Link object containing \code{linkfun}, \code{linkinv}, and
#' \code{mu.eta}.
#' @param spt Numeric vector of length two giving the lower and upper bounds of
#' the response support.
#' @param mu0 Numeric target mean used to identify the reference distribution.
#' @param init List of initial values for the MCMC algorithm.
#' @param cdpglmControl Object of class \code{"cdpglmControl"} containing MCMC,
#' tuning, and copula-control parameters.
#' @param thetaControl Object of class \code{"thetaControl"} containing control
#' arguments for the theta update procedure.
#'
#' @return
#' A list containing posterior samples, acceptance rates, prior settings,
#' support information, copula settings, and model specification information
#' used by \code{\link{cdpglm}}.
#'
#' @keywords internal
cdpglmFit <- function(formula, data, X, y, group_index,
                      link,
                      spt, mu0, init,
                      cdpglmControl, thetaControl) {
  # Extract cdpglmControl parameters
  burnin <- cdpglmControl$burnin
  thin <- cdpglmControl$thin
  save <- cdpglmControl$save
  mb <- cdpglmControl$mb
  Sb <- cdpglmControl$Sb
  M <- cdpglmControl$M
  alpha <- cdpglmControl$alpha
  delta <- cdpglmControl$delta
  c0 <- cdpglmControl$c0
  gamma <- cdpglmControl$gamma
  eps <- cdpglmControl$eps
  n_grid <- cdpglmControl$n_grid
  r_grid <- cdpglmControl$r_grid
  seed <- cdpglmControl$seed
  robust <- cdpglmControl$robust
  corr <- cdpglmControl$corr
  rho_proposal_sd <- cdpglmControl$rho_proposal_sd
  rho_prior_shape <- cdpglmControl$rho_prior_shape

  # Support
  if (is.unsorted(spt)) spt <- sort(spt)
  min_y <- spt[1]
  max_y <- spt[2]

  # MCMC Initialization
  n <- length(y)
  X <- as.matrix(X, nrow = n)
  p <- ncol(X)
  iter <- burnin + thin * save

  if (!is.null(seed)) set.seed(seed)
  if (is.null(c0)) c0 <- c0_silverman(y) / 4

  linkfun <- link$linkfun
  linkinv <- link$linkinv

  beta_samples <- matrix(NA, nrow = save, ncol = p)
  beta_names <- colnames(X)
  if (is.null(beta_names)) beta_names <- paste0("beta_", seq_len(p) - 1L)
  beta_names[1] <- "Intercept"
  colnames(beta_samples) <- beta_names
  theta_samples <- matrix(NA, nrow = save, ncol = n)
  h_samples <- matrix(NA, nrow = save, ncol = n)
  u_samples <- z_samples <- matrix(NA, nrow = save, ncol = n)
  rho_samples <- numeric(save)
  crm_samples <- vector("list", save)
  lnlik_samples <- numeric(save)

  beta <- as.numeric(init$beta)
  beta_rw_cov <- as.matrix(init$varbeta)
  if (!all(dim(beta_rw_cov) == c(p, p))) {
    stop("The beta random-walk covariance must be a p by p matrix.")
  }
  beta_rw_chol <- tryCatch(
    chol(beta_rw_cov),
    error = function(e) {
      stop("The beta random-walk covariance must be positive definite.")
    }
  )

  meanY_x <- as.numeric(linkinv(X %*% beta))

  z.tld <- z_tld <- as.numeric(init$crm$z.tld)
  J.tld <- J_tld <- as.numeric(init$crm$J.tld)

  theta <- theta_ref <- as.numeric(init$theta)
  bounds <- rho_bounds(corr, group_index)
  rho <- check_rho_start(init$rho, bounds, corr)
  rho_width <- diff(bounds)

  h <- log_copula_contribution_by_obs(
    y = y,
    group_index = group_index,
    rho = rho,
    corr = corr,
    crm.atoms = z.tld,
    crm.jumps = J.tld,
    theta = theta,
    c0 = c0,
    min_y = min_y,
    max_y = max_y
  )


  ord <- order(J_tld)[1:M]
  RL <- z_tld[ord]
  RJ <- J_tld[ord]

  z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta, min_y, max_y, eps)

  btheta <- b_theta(theta, z.tld, J.tld)
  T_vec <- exp(btheta)
  u <- rgamma(n, shape = 1, rate = T_vec)

  resampled_z <- resample_zstar(z)
  zstar <- resampled_z$zstar
  nstar <- resampled_z$nstar

  count_beta <- count_crm <- count_rho <- 0
  burning <- TRUE

  theta0 <- getTheta(
    spt = z.tld,
    f0 = J.tld,
    mu = mu0,
    sampprobs = NULL,
    ySptIndex = NULL,
    thetaStart = NULL
  )$theta

  if (any(theta0 > 50)) {
    theta0[theta0 > 50] <- 50
    warning("Capped some values of theta0 at 50. Consider increasing M.")
  }

  temp <- exp(theta0 * z.tld - max(theta0 * z.tld))
  Jtilt <- temp * J.tld
  Jtld_0 <- Jtilt / sum(Jtilt)

  if (robust) {
    J.tld <- Jtilt
    J.tld_ll <- Jtld_0
  } else {
    J.tld_ll <- Jtld_0
  }

  lnlik <- loglik(
    linkinv = linkinv,
    z = z,
    X = X,
    beta = beta,
    atoms = z.tld,
    jumps = J.tld_ll
  )

  # Auxiliary CRM
  z.tld_aug <- z.tld
  J.tld_aug <- J.tld

  # Beta prior
  if (is.null(mb)) {
    mprime <- spt[1] + (spt[2] - spt[1]) * 0.25
    Mprime <- spt[1] + (spt[2] - spt[1]) * 0.75
    gmprime <- linkfun(mprime)
    gMprime <- linkfun(Mprime)
    mid_spt <- (gmprime + gMprime) / 2
    mb <- c(mid_spt, rep(0, p - 1))
  } else if (length(mb) != p) {
    stop("length(mb) must match the number of covariates.")
  }

  if (is.null(Sb)) {
    mprime <- spt[1] + (spt[2] - spt[1]) * 0.25
    Mprime <- spt[1] + (spt[2] - spt[1]) * 0.75
    gmprime <- linkfun(mprime)
    gMprime <- linkfun(Mprime)
    sdX <- c(apply(as.matrix(X[, -1], nrow = n), 2, sd))
    Sbvec <- (gMprime - gmprime)^2 * c(100, (gamma / (2 * sdX))^2)
    Sb <- sqrt(Sbvec)
  } else if (length(Sb) != p) {
    stop("length(Sb) must match the number of betas.")
  } else if (any(Sb <= 0)) {
    stop("Sb must be positive definite.")
  }

  # Start MCMC loop
  for (itr in seq_len(iter)) {
    if (burning && itr > burnin) burning <- FALSE

    # Beta update
    theta <- getTheta(
      spt = z.tld, f0 = J.tld, mu = meanY_x, thetaStart = theta,
      thetaControl = thetaControl
    )$theta
    beta_ <- beta + drop(t(beta_rw_chol) %*% rnorm(p))
    mean_z_ <- as.numeric(linkinv(X %*% beta_))

    if (
      all(is.finite(mean_z_)) &&
        min(mean_z_) >= min(z.tld) && max(mean_z_) <= max(z.tld)
    ) {
      h <- log_copula_contribution_by_obs(
        y, group_index, rho, corr, z.tld, J.tld, theta, c0, min_y, max_y
      )
      cr_logpost_beta <- logpost_beta(
        beta, linkinv, z, X, z.tld, J.tld, mb, Sb, h, theta = theta
      )
      theta_ <- getTheta(
        spt = z.tld, f0 = J.tld, mu = mean_z_, thetaStart = theta,
        thetaControl = thetaControl
      )$theta
      h_ <- log_copula_contribution_by_obs(
        y, group_index, rho, corr, z.tld, J.tld, theta_, c0, min_y, max_y
      )
      pr_logpost_beta <- logpost_beta(
        beta_, linkinv, z, X, z.tld, J.tld, mb, Sb, h_, theta = theta_
      )
      log_acc_prob <- pr_logpost_beta - cr_logpost_beta

      if (!is.na(log_acc_prob) && log(runif(1)) < log_acc_prob) {
        beta <- beta_
        meanY_x <- mean_z_
        theta <- theta_
        h <- h_
        if (!burning) count_beta <- count_beta + 1
      }
    }

    # Theta update
    theta <- getTheta(
      spt = z.tld,
      f0 = J.tld,
      mu = meanY_x,
      sampprobs = NULL,
      ySptIndex = NULL,
      thetaStart = theta,
      thetaControl = thetaControl
    )$theta

    # h update
    h <- log_copula_contribution_by_obs(
      y = y,
      group_index = group_index,
      rho = rho,
      corr = corr,
      crm.atoms = z.tld,
      crm.jumps = J.tld,
      theta = theta,
      c0 = c0,
      min_y = min_y,
      max_y = max_y
    )

    # u update
    u <- sampler_u(u, zstar, nstar, theta, alpha, delta, min_y, max_y, eps, h)

    # CRM update
    crm_star <- crm_sampler(
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
      n_grid = n_grid,
      r_grid = r_grid
    )

    z.tld_star <- c(crm_star$RL, crm_star$zstar)
    J.tld_star <- c(crm_star$RJ, crm_star$Jstar)

    if (min(meanY_x) >= min(z.tld_star) && max(meanY_x) <= max(z.tld_star)) {
      # MH step
      theta_star <- getTheta(
        spt = z.tld_star,
        f0 = J.tld_star,
        mu = meanY_x,
        sampprobs = NULL,
        ySptIndex = NULL,
        thetaStart = theta,
        thetaControl = thetaControl
      )$theta

      h_star <- log_copula_contribution_by_obs(
        y = y,
        group_index = group_index,
        rho = rho,
        corr = corr,
        crm.atoms = z.tld_star,
        crm.jumps = J.tld_star,
        theta = theta_star,
        c0 = c0,
        min_y = min_y,
        max_y = max_y
      )

      u_star <- sampler_u(
        u, zstar, nstar, theta_star, alpha, delta, min_y, max_y, eps, h_star
      )
      crm_star_aug <- crm_sampler(
        M, u_star, zstar, nstar, theta_star, alpha, min_y, max_y, eps, h_star,
        n_grid = n_grid, r_grid = r_grid
      )
      z.tld_star_aug <- c(crm_star_aug$RL, crm_star_aug$zstar)
      J.tld_star_aug <- c(crm_star_aug$RJ, crm_star_aug$Jstar)

      b_thstar_crmstar <- b_theta(theta_star, z.tld_star, J.tld_star)
      b_th_crm <- b_theta(theta, z.tld, J.tld)
      b_thstar_crm <- b_theta(theta_star, z.tld, J.tld)
      b_th_crmstar <- b_theta(theta, z.tld_star, J.tld_star)

      b_thref_crmstar2 <- b_theta(theta_ref, z.tld_star_aug, J.tld_star_aug)
      b_thstar_crmstar2 <- b_theta(theta_star, z.tld_star_aug, J.tld_star_aug)
      b_thref_crm2 <- b_theta(theta_ref, z.tld_aug, J.tld_aug)
      b_thstar_crm2 <- b_theta(theta_star, z.tld_aug, J.tld_aug)

      log_r <- sum(
        (theta_star - theta) * z +
          (h_star - h) -
          b_thstar_crmstar + b_th_crm -
          b_thstar_crm + b_th_crmstar
      ) +
        sum(
          -b_thref_crmstar2 + b_thstar_crmstar2 +
            b_thref_crm2 - b_thstar_crm2
        )

      if (log(runif(1)) < log_r) {
        RL <- crm_star$RL
        RJ <- crm_star$RJ
        zstar <- crm_star$zstar
        Jstar <- crm_star$Jstar
        z.tld <- c(RL, zstar)
        J.tld <- c(RJ, Jstar)
        z.tld_aug <- z.tld_star_aug
        J.tld_aug <- J.tld_star_aug
        theta <- theta_star
        h <- h_star
        if (!burning) count_crm <- count_crm + 1
      }
    }

    # z update
    z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta, min_y, max_y, eps)

    # zstar and nstar update
    resampled_z <- resample_zstar(z)
    zstar <- resampled_z$zstar
    nstar <- resampled_z$nstar

    theta0 <- getTheta(
      spt = z.tld,
      f0 = J.tld,
      mu = mu0,
      sampprobs = NULL,
      ySptIndex = NULL,
      thetaStart = NULL
    )$theta

    if (any(theta0 > 50)) {
      theta0[theta0 > 50] <- 50
      warning("Capped some values of theta0 at 50.")
    }

    # rho update
    rho_unit <- (rho - bounds[1]) / rho_width
    rho_prop_unit <- plogis(qlogis(rho_unit) + rnorm(1, 0, rho_proposal_sd))
    rho_prop <- bounds[1] + rho_width * rho_prop_unit
    # log posterior
    logpost_current <- log_post_rho(
      rho, h, rho_prior_shape[1], rho_prior_shape[2], bounds
    )

    h_prop <- if (rho_prop > bounds[1] && rho_prop < bounds[2]) {
      log_copula_contribution_by_obs(
        y = y,
        group_index = group_index,
        rho = rho_prop,
        corr = corr,
        crm.atoms = z.tld,
        crm.jumps = J.tld,
        theta = theta,
        c0 = c0,
        min_y = min_y,
        max_y = max_y
      )
    } else {
      rep(-Inf, n)
    }

    rho_pr1 <- rho_prior_shape[1]
    rho_pr2 <- rho_prior_shape[2]
    logpost_prop <- log_post_rho(rho_prop, h_prop, rho_pr1, rho_pr2, bounds)

    log_jacobian <- log(rho_prop_unit) + log1p(-rho_prop_unit) -
      log(rho_unit) - log1p(-rho_unit)

    log_rho <- logpost_prop - logpost_current + log_jacobian

    if (!is.na(log_rho) && log(runif(1)) < log_rho) {
      if (!burning) count_rho <- count_rho + 1
      rho <- rho_prop
      h <- h_prop
    }

    temp <- exp(theta0 * z.tld - max(theta0 * z.tld))
    Jtilt <- temp * J.tld
    Jtld_0 <- Jtilt / sum(Jtilt)

    if (robust) {
      J.tld <- Jtilt
      J.tld_ll <- Jtld_0
    } else {
      J.tld_ll <- Jtld_0
    }

    lnlik <- loglik(
      linkinv = linkinv,
      z = z,
      X = X,
      beta = beta,
      atoms = z.tld,
      jumps = J.tld_ll
    )

    # Storage
    if (itr > burnin && (itr - burnin) %% thin == 0) {
      j <- (itr - burnin) / thin

      beta_samples[j, ] <- beta
      theta_samples[j, ] <- theta
      z_samples[j, ] <- z
      u_samples[j, ] <- u
      h_samples[j, ] <- h
      rho_samples[j] <- rho
      crm_samples[[j]] <- list(z.tld = z.tld, J.tld = J.tld)
      lnlik_samples[j] <- lnlik
    }
  }

  crm_samples <- data.frame(
    z.tld = I(lapply(crm_samples, `[[`, "z.tld")),
    J.tld = I(lapply(crm_samples, `[[`, "J.tld"))
  )

  samples <- list(
    z = z_samples,
    beta = beta_samples,
    theta = theta_samples,
    h = h_samples,
    rho = rho_samples,
    crm = crm_samples
  )

  list(
    samples = samples,
    mb = mb,
    Sb = Sb,
    spt = spt,
    mu0 = mu0,
    beta_acceptance = count_beta / (iter - burnin),
    crm_acceptance = count_crm / (iter - burnin),
    rho_acceptance = count_rho / (iter - burnin)
  )
}

#' Control arguments for the \code{dpglm} algorithm.
#'
#' This function returns control arguments for the \code{dpglm} algorithm.
#'
#' @param burnin Number of burn-in MCMC iterations. Defaults to 100.
#' @param thin Factor by which to thin MCMC iterations. Defaults to 10.
#' @param save Number of MCMC samples to return. Defaults to 1000.
#' @param mb Prior mean for beta. Defaults to NULL.
#' @param Sb Prior scale / diagonal variance entries for beta. Defaults to NULL.
#' @param M Number of random CRM atoms proposed at each CRM update. Defaults to 50.
#' @param alpha CRM mass parameter. Defaults to 1.
#' @param delta Tuning parameter for the \code{u} update. Defaults to 2.
#' @param c0 Bandwidth parameter used in the latent \code{z} update. Defaults to NULL.
#' @param gamma Shrinkage parameter for the default prior variance on \code{beta}.
#' Defaults to 1. Ignored if \code{Sb} is specified.
#' @param mu0 Mean of the reference distribution. Defaults to \code{NULL}, in which case
#' \code{mean(y)} is used internally.
#' @param spt Theoretical support for the response in the DP-GLM. This should be a
#' length-2 vector giving the lower and upper bounds. Defaults to \code{NULL}, in which
#' case the support is inferred internally from the data as \code{c(min(y)-eps, max(y)+eps)}.
#' @param eps Padding used when constructing default support bounds. Defaults to 1e-6.
#' @param n_grid Number of grid points used for the CRM jump-size grid. Defaults to 1501.
#' @param r_grid Number of grid points used for the CRM response-support grid. Defaults to 1501.
#' @param mc_size Number of Monte Carlo samples used to approximate \code{N(v)}
#' @param betaStart Initial value for \code{beta}. Defaults to NULL.
#' @param varbetaStart Initial value for the proposal covariance of \code{beta}. Defaults to NULL.
#' @param thetaStart Initial value for \code{theta}. Defaults to NULL.
#' @param crmStart Initial value for the CRM state, given as a list with components
#' \code{z.tld} and \code{J.tld}. Defaults to NULL.
#' @param seed Random seed. Defaults to NULL.
#' @param robust Provides numerical stability. Defaults to \code{TRUE}. See details for more.
#'
#' @return Object of S3 class "dpglmControl".
#'
#' @details Setting \code{robust = TRUE} will tilt the CRM weights at each MCMC iteration
#' to have the desired mean \code{mu0}. If \code{robust = FALSE}, these weights can
#' vary more in magnitude and lead to downstream numerical instability.
#'
#' @return An object of class \code{"dpglmControl"}; a list of control
#' parameters for the DPGLM fitting function.
#'
#' @examples
#' \dontrun{
#' ctrl <- dpglm.control(burnin = 100,
#'                       thin = 2,
#'                       save = 500,
#'                       spt = c(0, 1),
#'                       seed = 123)
#' }
#'
#' @export
dpglm.control <- function(burnin = 100L, thin = 10L, save = 1000L,
                          mb = NULL, Sb = NULL,
                          M = 50L, alpha = 1, delta = 2, c0 = NULL,
                          gamma = 1, mu0 = NULL, spt = NULL, eps = 1e-6,
                          n_grid = 1501L, r_grid = 1501L, mc_size = 1000L,
                          betaStart = NULL, varbetaStart = NULL,
                          thetaStart = NULL, crmStart = NULL,
                          seed = NULL, robust=TRUE)
{
  if (burnin < 0 || floor(burnin) != burnin) {
    stop("Number of burn-in samples must be an integer >= 0")
  }
  if (thin < 1 || floor(thin) != thin) {
    stop("Thin must be an integer >= 1")
  }
  if (save < 1 || floor(save) != save) {
    stop("Number of saved iterations must be an integer >= 1")
  }

  if (!is.numeric(n_grid) ||
      length(n_grid) != 1L ||
      !is.finite(n_grid) ||
      n_grid < 2 ||
      n_grid != floor(n_grid)) {
    stop("n_grid must be an integer >= 2.")
  }

  if (!is.numeric(r_grid) ||
      length(r_grid) != 1L ||
      !is.finite(r_grid) ||
      r_grid < 2 ||
      r_grid != floor(r_grid)) {
    stop("r_grid must be an integer >= 2.")
  }

  if (
    !is.numeric(mc_size) ||
    length(mc_size) != 1L ||
    !is.finite(mc_size) ||
    mc_size < 1 ||
    mc_size != floor(mc_size)
  ) {
    stop("mc_size must be an integer >= 1.")
  }

  n_grid <- as.integer(n_grid)
  r_grid <- as.integer(r_grid)
  mc_size <- as.integer(mc_size)

  if (!is.null(spt)) {
    if (!is.numeric(spt) || length(spt) != 2) {
      stop("For dpglm, spt must be a numeric vector of length 2 giving the support bounds.")
    }
    if (spt[1] >= spt[2]) {
      stop("For dpglm, spt must satisfy spt[1] < spt[2].")
    }
  }

  ctrl <- list(
    burnin       = burnin,
    thin         = thin,
    save         = save,
    mb           = mb,
    Sb           = Sb,
    M            = M,
    alpha        = alpha,
    delta        = delta,
    c0           = c0,
    gamma        = gamma,
    mu0          = mu0,
    spt          = spt,
    eps          = eps,
    n_grid       = n_grid,
    r_grid       = r_grid,
    mc_size      = mc_size,
    betaStart    = betaStart,
    varbetaStart = varbetaStart,
    thetaStart   = thetaStart,
    crmStart     = crmStart,
    seed         = seed,
    robust       = robust
  )

  class(ctrl) <- "dpglmControl"
  ctrl
}


#' Main MCMC Fitting Function for DPGLM
#'
#' Performs MCMC to fit the DPGLM model.
#'
#' @param formula Model formula supplied to \code{\link{dpglm}}.
#' @param data Data frame containing the variables appearing in \code{formula}.
#' @param X Numeric design matrix.
#' @param y Numeric response vector.
#' @param link Link object containing \code{linkfun}, \code{linkinv}, and \code{mu.eta}.
#' @param spt Numeric vector of length two giving the lower and upper bounds
#' of the response support.
#' @param mu0 Numeric target mean used to identify the reference distribution.
#' @param init List of initial values for the MCMC algorithm.
#' @param dpglmControl Object of class \code{"dpglmControl"} containing MCMC
#' and tuning parameters.
#' @param thetaControl Object of class \code{"thetaControl"} containing control
#' arguments for the theta update procedure.
#'
#' @return
#' A list containing posterior samples, acceptance rates, prior settings,
#' support information, and model specification information used by
#' \code{\link{dpglm}}.
#'
#' @keywords internal
dpglmFit <- function(formula, data, X, y,        # Data
                     link,                       # Link
                     spt, mu0, init,             # Model Specs
                     dpglmControl, thetaControl) # Controls
{
  timer_now <- function() proc.time()[["elapsed"]]

  chunk_times <- c(
    beta_mh = 0,
    u_update = 0,
    zstar_resample = 0,
    crm_sampler = 0,
    crm_accept = 0,
    z_update = 0,
    storage = 0
  )

  # Extract dpglmControl parameters
  burnin       <- dpglmControl$burnin
  thin         <- dpglmControl$thin
  save         <- dpglmControl$save
  mb           <- dpglmControl$mb
  Sb           <- dpglmControl$Sb
  M            <- dpglmControl$M
  alpha        <- dpglmControl$alpha
  delta        <- dpglmControl$delta
  c0           <- dpglmControl$c0
  gamma        <- dpglmControl$gamma
  eps          <- dpglmControl$eps
  n_grid       <- dpglmControl$n_grid
  r_grid       <- dpglmControl$r_grid
  mc_size      <- dpglmControl$mc_size
  seed         <- dpglmControl$seed
  robust       <- dpglmControl$robust

  # Set theoretical support bounds
  if (is.unsorted(spt)) spt <- sort(spt)
  min_y <- spt[1]
  max_y <- spt[2]

  # MCMC Initialization
  n    <- length(y)
  X    <- as.matrix(X, nrow=n)
  p    <- dim(X)[2]
  iter <- burnin + thin * save

  # Default h for the DPGLM model
  h <- h_ <- hstar <- 0

  if (!is.null(seed)) set.seed(seed)

  if (is.null(c0)) c0 <- c0_silverman(y) / 4

  # Extract link
  linkfun <- link$linkfun
  linkinv <- link$linkinv

  beta_samples <- matrix(NA, nrow = save, ncol = p)
  beta_names <- colnames(X)
  if (is.null(beta_names)) beta_names <- paste0("beta_", seq_len(p) - 1L)
  beta_names[1] <- "Intercept"
  colnames(beta_samples) <- beta_names
  theta_samples <- matrix(NA, nrow = save, ncol = n)
  u_samples <- z_samples <- matrix(NA, nrow = save, ncol = n)
  crm_samples   <- vector("list", save)
  ######
  crm_u_samples <- matrix(NA_real_, nrow = save, ncol = n)
  crm_theta_samples <- matrix(NA_real_, nrow = save, ncol = n)
  #####
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
  meanY_x <- linkinv(X %*% beta)
  z.tld <- z_tld <-  as.numeric(init$crm$z.tld)
  J.tld <- J_tld <- as.numeric(init$crm$J.tld)

  ord <- order(J_tld)[1:M]

  RL <- z_tld[ord]
  RJ <- J_tld[ord]
  theta <- theta_ref <- init$theta
  z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta, min_y, max_y, eps)

  btheta <- b_theta(theta, z.tld, J.tld)

  T_vec <- exp(btheta)
  u <- rgamma(n, shape = 1, rate = T_vec)

  resampled_z <- resample_zstar(z)
  zstar <- resampled_z$zstar
  nstar <- resampled_z$nstar
  Jstar <- rgamma(n = length(nstar), shape = nstar, rate = 1)

  count1 <- count2 <- 0
  burning <- TRUE

  theta0 <- getTheta(
    spt = z.tld,
    f0  = J.tld,
    mu  = mu0,
    sampprobs  = NULL,
    ySptIndex  = NULL,
    thetaStart = NULL
  )$theta

  if (any(theta0 > 50)) {
    idx <- which(theta0 > 50)
    theta0[idx] <- 50
    warning("Capped some values of Theta at 50. Consider increasing M.")
  }

  temp <- exp(theta0 * z.tld - max(theta0 * z.tld))
  Jtilt <- temp * J.tld
  W <- sum(Jtilt)
  Jtld_0 <- Jtilt / W

  if (robust) {
    J.tld <- Jtilt
    J.tld_ll <- Jtld_0
  } else {
    J.tld_ll <- Jtld_0
  }

  lnlik <- loglik(linkinv = linkinv, z = z, X = X, beta = beta, atoms = z.tld, jumps = J.tld_ll)

  # Beta prior
  if (is.null(mb)) {
    mprime  <- spt[1] + (spt[2] - spt[1]) * 0.25
    Mprime  <- spt[1] + (spt[2] - spt[1]) * 0.75
    gmprime <- linkfun(mprime)
    gMprime <- linkfun(Mprime)
    mid_spt <- (gmprime + gMprime) / 2
    mb <- c(mid_spt, rep(0, p-1))
  }
  else if (length(mb) != p) stop("length(mb) must match the number of covariates.")

  if (is.null(Sb)) {
    mprime  <- spt[1] + (spt[2] - spt[1]) * 0.25
    Mprime  <- spt[1] + (spt[2] - spt[1]) * 0.75
    gmprime <- linkfun(mprime)
    gMprime <- linkfun(Mprime)
    sdX     <- c(apply(as.matrix(X[, -1], nrow=n), 2, sd))
    Sbvec   <- (gMprime - gmprime)^2 * c(100, (gamma / (2 * sdX))^2) # Re-scaling on the linear-predictor scale
    Sb <- sqrt(Sbvec) # Entries corresponding to sqrt(diagonal) of beta prior cov
  } else if (length(Sb) != p) stop("length(Sb) must match the number of betas.")
  else if (any(Sb <= 0)) stop("Sb must be positive definite.")

  # Start MCMC loop
  for(itr in seq_len(iter)){
    if (burning && itr > burnin) burning <- FALSE

    #######
    save_now <- (itr > burnin &&(itr - burnin) %% thin == 0L)
    if (save_now) save_index <- as.integer((itr - burnin) / thin)
    #######


    # beta update ------------------------------------
    t0 <- timer_now()
    beta_ <- beta + drop(t(beta_rw_chol) %*% rnorm(p))
    mean_z_ <- as.numeric(linkinv(X %*% beta_))

    proposal_valid <- (
      all(is.finite(mean_z_)) &&
        min(mean_z_) >= min(z.tld) &&
        max(mean_z_) <= max(z.tld)
    )

    if (proposal_valid) {

      cr_logpost_beta <- logpost_beta(
        beta = beta,
        linkinv = linkinv,
        z = z,
        X = X,
        atoms = z.tld,
        jumps = J.tld,
        mu_beta = mb,
        sigma_beta = Sb,
        h = h
      )

      pr_logpost_beta <- logpost_beta(
        beta = beta_,
        linkinv = linkinv,
        z = z,
        X = X,
        atoms = z.tld,
        jumps = J.tld,
        mu_beta = mb,
        sigma_beta = Sb,
        h = h
      )

      log_acc_prob <- pr_logpost_beta - cr_logpost_beta

      if (
        is.finite(log_acc_prob) &&
        log(runif(1)) < min(0, log_acc_prob)
      ) {
        beta <- beta_
        meanY_x <- mean_z <- mean_z_

        if (!burning) {
          count1 <- count1 + 1
        }
      }
    }

    chunk_times["beta_mh"] <- chunk_times["beta_mh"] + timer_now() - t0

    # theta update ------------------------------------
    theta <- getTheta(
      spt = z.tld,
      f0  = J.tld,
      mu  = meanY_x,
      sampprobs  = NULL,
      ySptIndex  = NULL,
      thetaStart = theta
    )$theta

    if (any(theta > 50)) {
      idx <- which(theta > 50)
      theta[idx] <- 50
      warning("Capped some values of Theta at 50. Consider increasing M.")
    }

    # u update ----------------------------------------
    t0 <- timer_now()
    u <- sampler_u(u, zstar, nstar, theta, alpha, delta, min_y, max_y, eps, h)
    chunk_times["u_update"] <- chunk_times["u_update"] + timer_now() - t0

    # CRM update --------------------------------------
    t0 <- timer_now()
    ########
    if (save_now) {
      crm_u_samples[save_index, ] <- u
      crm_theta_samples[save_index, ] <- theta
    }
    ########
    crm_draws <- crm_sampler(
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
      n_draws = 2L,
      n_grid = n_grid,
      r_grid = r_grid,
      mc_size = mc_size
    )

    crm_star <- crm_draws[[1]]
    crm_2    <- crm_draws[[2]]

    z.tld_star <- c(crm_star$RL, crm_star$zstar)
    J.tld_star <- c(crm_star$RJ, crm_star$Jstar)
    z.tld_2 <- c(crm_2$RL, crm_2$zstar)
    J.tld_2 <- c(crm_2$RJ, crm_2$Jstar)
    chunk_times["crm_sampler"] <- chunk_times["crm_sampler"] + timer_now() - t0

    if(min(meanY_x) >= min(z.tld_star) && max(meanY_x) <= max(z.tld_star)){
      # MH step
      theta_star <- getTheta(
        spt = z.tld_star,
        f0  = J.tld_star,
        mu  = meanY_x,
        sampprobs  = NULL,
        ySptIndex  = NULL,
        thetaStart = theta
      )$theta

      if (any(theta_star > 50)) {
        idx <- which(theta_star > 50)
        theta_star[idx] <- 50
        warning("Capped some values of Theta at 50. Consider increasing M.")
      }

      t0 <- timer_now()
      crm_star_2 <- crm_sampler(
        M = M,
        u = u,
        zstar = zstar,
        nstar = nstar,
        theta = theta_star,
        alpha = alpha,
        min_y = min_y,
        max_y = max_y,
        eps = eps,
        h = h,
        n_draws = 1L,
        n_grid = n_grid,
        r_grid = r_grid,
        mc_size = mc_size
      )
      z.tld_star_2 <- c(crm_star_2$RL, crm_star_2$zstar)
      J.tld_star_2 <- c(crm_star_2$RJ, crm_star_2$Jstar)
      chunk_times["crm_sampler"] <- chunk_times["crm_sampler"] + timer_now() - t0

      t0 <- timer_now()
      b_thstar_crmstar <- b_theta(theta_star, z.tld_star, J.tld_star)
      b_th_crm         <- b_theta(theta,      z.tld,      J.tld)
      b_thstar_crm     <- b_theta(theta_star, z.tld,      J.tld)
      b_th_crmstar     <- b_theta(theta,      z.tld_star, J.tld_star)

      b_th_crmstar2    <- b_theta(theta_ref,      z.tld_star_2, J.tld_star_2)
      b_thstar_crmstar2 <- b_theta(theta_star, z.tld_star_2, J.tld_star_2)
      b_th_crm2        <- b_theta(theta_ref,      z.tld_2,      J.tld_2)
      b_thstar_crm2    <- b_theta(theta_star, z.tld_2,      J.tld_2)

      log_r <- sum(
        (theta_star - theta) * z -
          b_thstar_crmstar + b_th_crm -
          b_thstar_crm + b_th_crmstar
      ) +
        sum(
          -b_th_crmstar2 + b_thstar_crmstar2 +
            b_th_crm2 - b_thstar_crm2
        )

      if(log(runif(1)) < log_r){
        RL <- crm_star$RL
        RJ <- crm_star$RJ
        zstar <- crm_star$zstar
        Jstar <- crm_star$Jstar
        z.tld <- c(RL, zstar)
        J.tld <- c(RJ, Jstar)
        theta <- theta_star
        if (!burning) count2 <- count2 + 1
      }
      chunk_times["crm_accept"] <- chunk_times["crm_accept"] + timer_now() - t0
    }

    # z update ------------------------------------------------------------------
    t0 <- timer_now()
    z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta, min_y, max_y, eps)
    chunk_times["z_update"] <- chunk_times["z_update"] + timer_now() - t0

    # zstar and nstar update ----------------------------------------------------
    t0 <- timer_now()
    resampled_z <- resample_zstar(z)
    chunk_times["zstar_resample"] <- chunk_times["zstar_resample"] + timer_now() - t0
    zstar <- resampled_z$zstar
    nstar <- resampled_z$nstar

    theta0 <- getTheta(
      spt = z.tld,
      f0  = J.tld,
      mu  = mu0,
      sampprobs  = NULL,
      ySptIndex  = NULL,
      thetaStart = NULL
    )$theta

    if (any(theta0 > 50)) {
      idx <- which(theta0 > 50)
      theta0[idx] <- 50
      warning("Capped some values of Theta at 50. Consider increasing M.")
    }

    temp <- exp(theta0 * z.tld - max(theta0 * z.tld))
    Jtilt <- temp * J.tld
    W <- sum(Jtilt)
    Jtld_0 <- Jtilt / W
    if (robust) {
    J.tld <- Jtilt
    J.tld_ll <- Jtld_0
    } else {
      J.tld_ll <- Jtld_0
    }

    lnlik <- loglik(linkinv = linkinv, z = z, X = X, beta = beta, atoms = z.tld, jumps = J.tld_ll)

    # Storing MCMC simulations --------------------------------------------------
    t0 <- timer_now()
    if (itr > burnin && (itr - burnin) %% thin == 0) {
      j <- (itr - burnin) / thin

      beta_samples[j, ] <- beta
      theta_samples[j, ] <- theta
      z_samples[j, ] <- z
      u_samples[j, ] <- u
      crm_samples[[j]] <- list(z.tld = z.tld, J.tld = J.tld)
      lnlik_samples[j] <- lnlik
    }
    chunk_times["storage"] <- chunk_times["storage"] + timer_now() - t0
  }

  crm_samples <- data.frame(
  z.tld = I(lapply(crm_samples, `[[`, "z.tld")),
  J.tld = I(lapply(crm_samples, `[[`, "J.tld"))
  )

  samples <- list(z = z_samples, beta = beta_samples,
                  crm = crm_samples, theta = theta_samples)

  #########
  crm_inputs <- list(
    u = crm_u_samples,
    theta = crm_theta_samples
  )
  #########

  timing <- data.frame(
    chunk = names(chunk_times),
    seconds = as.numeric(chunk_times),
    proportion = as.numeric(chunk_times) / sum(chunk_times)
  )

  list(samples         = samples,
       crm_inputs      = crm_inputs,
       mb              = mb,
       Sb              = Sb,
       spt             = spt,
       mu0             = mu0,
       beta_acceptance = count1 / (iter - burnin),
       crm_acceptance  = count2 / (iter - burnin),
       timing          = timing)
}

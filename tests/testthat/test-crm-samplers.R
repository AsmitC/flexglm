test_that("CRM draws and RNG consumption match the peer reference", {
  ref <- new.env()
  sys.source(test_path("fixtures", "peer-samplers.R"), ref)
  for (case in ref$peer_sampler_fixtures$crm[c(
    "thinning", "inversion", "multiple_draws"
  )]) {
    set.seed(case$seed)
    actual <- do.call(crm_sampler, case$args)
    expect_equal(actual, case$expected, tolerance = 1e-12)
    expect_equal(runif(3), case$next_random, tolerance = 1e-14)
  }
})

test_that("ten-sweep auxiliary updates match the peer reference", {
  ref <- new.env()
  sys.source(test_path("fixtures", "peer-samplers.R"), ref)
  for (case in ref$peer_sampler_fixtures$u) {
    # Both C++ samplers update the supplied u vector in place.
    args <- unserialize(serialize(case$args, NULL))
    set.seed(case$seed)
    actual <- do.call(sampler_u, args)
    expect_equal(actual, case$expected, tolerance = 1e-12)
    expect_equal(runif(3), case$next_random, tolerance = 1e-14)
    expect_equal(args$u, case$input_after, tolerance = 1e-12)
    expect_equal(do.call(log_post_u, case$args[names(case$args) != "delta"]),
      case$log_post,
      tolerance = 1e-12
    )
  }
})

test_that("CRM grid caching respects controls and scalar h is supported", {
  args <- list(
    M = 10L, u = c(1, 2), zstar = c(0.3, 0.7), nstar = c(1, 1),
    theta = c(-1, 1), alpha = 1, min_y = 0, max_y = 1,
    eps = 1e-6, h = 0, n_grid = 71L, r_grid = 53L
  )
  set.seed(71)
  original <- do.call(crm_sampler, args)
  changed <- args
  changed$n_grid <- 101L
  changed$eps <- 2e-6
  other <- do.call(crm_sampler, changed)
  grid <- -log(seq(exp(-changed$eps), exp(-5e-4), length.out = changed$n_grid))
  expect_true(all(vapply(
    other$RJ, function(x) min(abs(grid - x)) < 1e-14,
    logical(1)
  )))
  args$h <- c(0, 0)
  set.seed(71)
  expect_identical(do.call(crm_sampler, args), original)
  expect_true(all(other$RL >= changed$min_y + changed$eps))
  expect_true(all(other$RL <= changed$max_y - changed$eps))

  u_args <- args[names(args) %in% c(
    "u", "zstar", "nstar", "theta", "alpha",
    "min_y", "max_y", "eps", "h"
  )]
  u_args$delta <- 2
  set.seed(71)
  vector_h <- do.call(sampler_u, unserialize(serialize(u_args, NULL)))
  u_args$h <- 0
  set.seed(71)
  expect_identical(do.call(sampler_u, u_args), vector_h)
})

test_that("CRM grid controls agree across DPGLM and CDPGLM", {
  for (control in list(dpglm.control, cdpglm.control)) {
    ctrl <- control()
    expect_equal(ctrl$n_grid, 3001L)
    expect_equal(ctrl$r_grid, 3001L)
    expect_equal(control(n_grid = 71, r_grid = 53)$n_grid, 71L)
    expect_error(control(n_grid = 1), "n_grid")
    expect_error(control(r_grid = 2.5), "r_grid")
    expect_false("mc_size" %in% names(ctrl))
    expect_false("mc_size" %in% names(formals(control)))
  }
  expect_false(any(c("n_grid", "r_grid") %in% names(dirglm.control())))
  expect_false("mc_size" %in% names(formals(crm_sampler)))
  expect_false("mc_size" %in% names(formals(crm_sampler_cpp)))
})

# Rescale alpha to compare support-width quadrature with the peer reference.
test_that("CRM quadrature includes the response support width", {
  ref <- new.env()
  sys.source(test_path("fixtures", "peer-samplers.R"), ref)
  case <- ref$peer_sampler_fixtures$crm$shifted_support
  args <- case$args
  args$alpha <- args$alpha * (1 - 2 * args$eps) /
    (args$max_y - args$min_y - 2 * args$eps)
  set.seed(case$seed)
  expect_equal(do.call(crm_sampler, args), case$expected, tolerance = 1e-12)
  expect_equal(runif(3), case$next_random, tolerance = 1e-14)
})

test_that("both fitting drivers retain the auxiliary CRM across rejection", {
  for (model in c("dpglm", "cdpglm")) {
    fitter <- get(paste0(model, "Fit"))
    env <- new.env(parent = environment(fitter))
    environment(fitter) <- env
    trace <- new.env(parent = emptyenv())
    trace$draws <- list()
    trace$normalizers <- list()
    trace$u_calls <- list()
    trace$b_index <- 0L

    # Distinct current/proposed values expose accidentally reused theta or h.
    env$getTheta <- function(spt, f0, mu, ...) {
      list(theta = rep(f0[1] / 100, length(mu)))
    }
    env$loglik <- function(...) 0
    env$logpost_beta <- function(...) 0
    env$optim <- function(...) NULL
    env$log_copula_contribution_by_obs <- function(y, crm.jumps, ...) {
      rep(sum(crm.jumps) / 1000, length(y))
    }
    env$z_sampler_unifK <- function(y, ...) rep(mean(y), length(y))
    env$sampler_u <- function(u, zstar, nstar, theta, alpha, delta,
                              min_y, max_y, eps, h) {
      trace$u_calls[[length(trace$u_calls) + 1L]] <- list(u = u,
                                                          theta = theta, h = h)
      u + 1
    }
    env$crm_sampler <- function(M, u, zstar, nstar, theta, alpha,
                                min_y, max_y, eps, h, n_draws = 1L, ...) {
      i <- length(trace$draws) + 1L
      iteration <- (i + 1L) %/% 2L
      weight <- if (i %% 2L) 10 * iteration else 100 * iteration
      draw <- list(
        RL = c(0.1, 0.9), RJ = c(weight, weight + 1),
        zstar = numeric(), Jstar = numeric()
      )
      trace$draws[[i]] <- list(
        u = u, theta = theta, h = h,
        n_draws = n_draws, draw = draw
      )
      trace$b_index <- 0L
      draw
    }
    env$b_theta <- function(theta, spt, f0) {
      if (!length(trace$draws)) {
        return(rep(0, length(theta)))
      }
      trace$b_index <- trace$b_index + 1L
      iteration <- length(trace$draws) %/% 2L
      trace$normalizers[[length(trace$normalizers) + 1L]] <-
        list(theta = theta, spt = spt, f0 = f0)
      # Force the second CRM proposal to be rejected; accept the first and
      # third.
      rep(if (iteration == 2L && trace$b_index == 1L) 1e9 else 0, length(theta))
    }

    y <- c(0.45, 0.48, 0.52, 0.55)
    X <- cbind(Intercept = 1, x = c(-1.5, -0.5, 0.5, 1.5))
    init <- list(
      beta = c(0, 0), varbeta = diag(0.01, 2), theta = rep(0, 4),
      crm = list(z.tld = c(0.1, 0.9), J.tld = c(0.5, 0.5)), rho = 0.3
    )
    ctrl <- do.call(
      get(paste0(model, ".control")),
      list(
        burnin = 0, thin = 1, save = 3, M = 2, c0 = 0.2,
        mb = c(0, 0), Sb = c(1, 1), robust = FALSE, seed = 71
      )
    )
    args <- list(
      formula = y ~ x, data = data.frame(y, x = X[, 2]),
      X = X, y = y, link = make.link("logit"), spt = c(0, 1),
      mu0 = 0.5, init = init, thetaControl = theta.control()
    )
    args[[paste0(model, "Control")]] <- ctrl
    if (model == "cdpglm") args$group_index <- rep(1:2, each = 2)
    fit <- do.call(fitter, args)

    expect_length(trace$draws, 6L)
    expect_length(trace$u_calls, 6L)
    expect_true(all(vapply(trace$draws, function(x) x$n_draws == 1L,
                           logical(1))))
    for (i in c(1L, 3L, 5L)) {
      expect_equal(trace$draws[[i + 1L]]$u, trace$draws[[i]]$u + 1)
      for (j in c(i, i + 1L)) {
        expect_equal(trace$draws[[j]]$theta, trace$u_calls[[j]]$theta)
        expect_equal(trace$draws[[j]]$h, trace$u_calls[[j]]$h)
      }
    }
    expect_length(trace$normalizers, 24L)
    for (iteration in seq_len(3L)) {
      current <- if (iteration == 1L) c(0.5, 0.5) else c(10, 11)
      proposed <- c(10 * iteration, 10 * iteration + 1)
      retained_aux <- if (iteration == 1L) c(0.5, 0.5) else c(100, 101)
      proposed_aux <- c(100 * iteration, 100 * iteration + 1)
      current_theta <- rep(current[1] / 100, length(y))
      proposed_theta <- rep(proposed[1] / 100, length(y))
      i <- 2L * iteration - 1L
      expect_equal(trace$draws[[i]]$theta, current_theta)
      expect_equal(trace$draws[[i + 1L]]$theta, proposed_theta)
      expected_h <- function(weights) {
        if (model == "cdpglm") rep(sum(weights) / 1000, length(y)) else 0
      }
      expect_equal(trace$draws[[i]]$h, expected_h(current))
      expect_equal(trace$draws[[i + 1L]]$h, expected_h(proposed))

      # All eight normalizers must use the peer's current/proposed/reference
      # theta and main/auxiliary CRM pairings, including the retained state.
      expected_weights <- list(
        proposed, current, current, proposed,
        proposed_aux, proposed_aux, retained_aux, retained_aux
      )
      expected_theta <- list(
        proposed_theta, current_theta, proposed_theta, current_theta,
        init$theta, proposed_theta, init$theta, proposed_theta
      )
      for (j in seq_len(8L)) {
        actual <- trace$normalizers[[(iteration - 1L) * 8L + j]]
        expect_equal(actual$f0, expected_weights[[j]])
        expect_equal(actual$theta, expected_theta[[j]])
        expect_equal(actual$spt, c(0.1, 0.9))
      }
    }
    expect_equal(fit$samples$crm$J.tld, I(list(c(10, 11), c(10, 11), c(30,
                                                                       31))))
    expect_equal(fit$crm_acceptance, 2 / 3)
  }
})

test_that("robust stabilization and fitted interfaces survive the migration", {
  dat <- data.frame(
    y = c(0.42, 0.45, 0.48, 0.50, 0.52, 0.55, 0.58, 0.60),
    x = c(-1.5, -1, -0.5, 0, 0.25, 0.5, 1, 1.5)
  )
  atoms <- seq(0.01, 0.99, length.out = 25)
  common <- list(
    burnin = 1, thin = 1, save = 3, M = 5, spt = c(0, 1),
    mu0 = mean(dat$y), c0 = 0.15, mb = c(0, 0), Sb = c(1, 1),
    betaStart = c(0, 0), varbetaStart = diag(0.01, 2),
    thetaStart = rep(0, nrow(dat)),
    crmStart = list(z.tld = atoms, J.tld = rep(1 / 25, 25)), seed = 832
  )

  for (model in c("dpglm", "cdpglm")) {
    for (robust in c(TRUE, FALSE)) {
      controls <- c(common, list(robust = robust))
      if (model == "cdpglm") {
        controls$corr <- "ar1"
        controls$rhoStart <- 0.3
      }
      ctrl <- do.call(get(paste0(model, ".control")), controls)
      args <- list(formula = y ~ x, data = dat, link = "logit")
      args[[paste0(model, "Control")]] <- ctrl
      if (model == "cdpglm") args$group_index <- rep(1:4, each = 2)
      fit <- suppressWarnings(do.call(get(model), args))
      repeated <- suppressWarnings(do.call(get(model), args))
      expect_equal(fit$samples, repeated$samples)
      expect_s3_class(fit, model)
      expect_true(all(is.finite(fit$samples$beta)))
      expect_true(all(is.finite(fit$samples$theta)))
      if (robust) {
        means <- mapply(
          function(z, J) sum(z * J) / sum(J),
          fit$samples$crm$z.tld, fit$samples$crm$J.tld
        )
        expect_equal(as.numeric(means), rep(common$mu0, common$save),
                     tolerance = 1e-7)
      }
      if (model == "dpglm") {
        expect_equal(dim(fit$crm_inputs$u), c(3, nrow(dat)))
        expect_equal(dim(fit$crm_inputs$theta), c(3, nrow(dat)))
        expect_true(all(is.finite(fit$crm_inputs$u)))
        expect_true(all(fit$timing$seconds >= 0))
      }
    }
  }
})

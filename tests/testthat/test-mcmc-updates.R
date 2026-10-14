update_test_data <- function() {
  data.frame(y = c(0.42, 0.45, 0.48, 0.50, 0.52, 0.55, 0.58, 0.60),
             x = c(-1.5, -1, -0.5, 0, 0.25, 0.5, 1, 1.5),
             id = rep(1:4, each = 2))
}

update_test_control <- function() {
  atoms <- seq(0.01, 0.99, length.out = 25)
  list(burnin = 0, thin = 1, save = 2, spt = c(0, 1), c0 = 0.15,
       M = 5, betaStart = c(0, 0), varbetaStart = diag(0.01, 2),
       thetaStart = rep(0, 8), mb = c(0, 0), Sb = c(1, 1),
       crmStart = list(z.tld = atoms, J.tld = rep(1 / 25, 25)), seed = 832)
}

test_that("rho bounds follow the correlation structure and group sizes", {
  groups <- c(rep(1, 4), rep(2, 2))
  expect_equal(rho_bounds("ar1", groups), c(-1, 1))
  expect_equal(rho_bounds("ex", groups), c(-1 / 3, 1))
  expect_equal(rho_bounds("ex", seq_len(4)), c(-1, 1))
  expect_equal(check_rho_start(NULL, c(-1, 1), "ar1"), 0)
  expect_equal(check_rho_start(NULL, c(-1 / 3, 1), "ex"), 1 / 3)
  expect_warning(value <- check_rho_start(-0.5, c(-1 / 3, 1), "ex"), "rhoStart")
  expect_equal(value, 1 / 3)
  expect_equal(check_rho_start(-0.5, c(-1, 1), "ar1"), -0.5)
  expect_warning(ctrl <- cdpglm.control(rhoStart = 1), "rhoStart")
  expect_null(ctrl$rhoStart)

  bounds <- c(-1 / 3, 1)
  unit_rho <- 0.25
  rho <- bounds[1] + diff(bounds) * unit_rho
  expect_equal(log_post_rho(rho, c(-1, -2), 8, 2, bounds),
               -3 + dbeta(unit_rho, 8, 2, log = TRUE) - log(diff(bounds)))
  expect_equal(log_post_rho(bounds[1], 0, 8, 2, bounds), -Inf)
  expect_equal(log_post_rho(bounds[2], 0, 8, 2, bounds), -Inf)
})

test_that("CDPGLM supports valid negative correlations", {
  dat <- update_test_data()
  for (corr in c("ex", "ar1")) {
    controls <- update_test_control()
    controls$corr <- corr
    controls$rhoStart <- -0.4
    controls$n_grid <- 71L
    controls$r_grid <- 53L
    fit <- suppressWarnings(cdpglm(
      y ~ x, dat, group_index = dat$id,
      cdpglmControl = do.call(cdpglm.control, controls)
    ))
    bounds <- rho_bounds(corr, dat$id)
    expect_true(all(is.finite(fit$samples$rho)))
    expect_true(all(fit$samples$rho > bounds[1]))
    expect_true(all(fit$samples$rho < bounds[2]))
    expect_true(any(fit$samples$rho < 0))
  }
})

test_that("beta random walks condition on the current imputed state", {
  dat <- update_test_data()
  X <- model.matrix(y ~ x, dat)
  for (model in c("dpglm", "cdpglm")) {
    fitter <- get(paste0(model, "Fit"))
    env <- new.env(parent = environment(fitter))
    environment(fitter) <- env
    calls <- list()
    env$optim <- function(...) stop("Beta updates must not optimize.")
    env$rnorm <- function(n, ...) {
      if (n == 2L) c(0.1, -0.2) else stats::rnorm(n, ...)
    }
    env$logpost_beta <- function(beta, linkinv, z, X, atoms, jumps,
                                 mu_beta, sigma_beta, h, theta = NULL) {
      calls[[length(calls) + 1L]] <<-
        list(beta = beta, z = z, atoms = atoms, jumps = jumps, h = h)
      length(calls)
    }
    atoms <- seq(0.01, 0.99, length.out = 25)
    jumps <- rep(1 / 25, 25)
    init <- list(beta = c(0, 0), varbeta = diag(0.01, 2),
                 theta = rep(0, nrow(dat)),
                 crm = list(z.tld = atoms, J.tld = jumps), rho = 0.3)
    ctrl <- do.call(get(paste0(model, ".control")), list(
      burnin = 0, thin = 1, save = 1, M = 5, c0 = 0.15,
      mb = c(0, 0), Sb = c(1, 1), n_grid = 71, r_grid = 53,
      robust = FALSE, seed = 832
    ))
    args <- list(formula = y ~ x, data = dat, X = X, y = dat$y,
                 link = make.link("logit"), spt = c(0, 1), mu0 = mean(dat$y),
                 init = init, thetaControl = theta.control())
    args[[paste0(model, "Control")]] <- ctrl
    if (model == "cdpglm") args$group_index <- dat$id
    fit <- do.call(fitter, args)
    steps <- if (model == "dpglm") 2L else 1L
    expect_length(calls, 2L * steps)
    expect_equal(as.numeric(fit$samples$beta[1, ]), steps * c(0.01, -0.02))
    for (i in seq_along(calls)) {
      call <- calls[[i]]
      expect_equal(call$atoms, atoms)
      expect_equal(call$jumps, jumps)
      expect_equal(call$z, calls[[1]]$z)
      if (model == "cdpglm") {
        theta <- getTheta(atoms, jumps, plogis(X %*% call$beta))$theta
        expected <- log_copula_contribution_by_obs(
          dat$y, dat$id, 0.3, "ex", atoms, jumps, theta, 0.15, 0, 1
        )
        expect_equal(call$h, expected, tolerance = 1e-7)
      } else {
        expect_equal(call$h, 0)
      }
    }
  }
})

test_that("componentwise DirGLM updates refresh scales and count acceptances", {
  updater <- beta_update_separate
  env <- new.env(parent = environment(updater))
  environment(updater) <- env
  scales <- numeric()
  draws <- 0L
  env$rnorm <- function(n, mean, sd) {
    scales <<- c(scales, sd)
    1
  }
  env$tht_sol <- function(...) {
    draws <<- draws + 1L
    list(tht = c(0, 0), bpr2 = c(1, 1),
         btht = rep(if (draws == 1L) -1000 else 1000, 2))
  }
  env$Sigma_beta <- function(...) diag(c(4, 9))
  out <- updater(
    diag(2), c(0.4, 0.6), c(0, 1), c(0, 0), diag(2), c(0.5, 0.5),
    c(0, 0), c(1, 1), c(0, 0), 0.1, qlogis, plogis,
    make.link("logit")$mu.eta, c(0, 0), c(100, 100)
  )
  expect_equal(scales, c(1, 3))
  expect_equal(out$cr_bt, c(1, 0))
  expect_equal(out$acc_beta, 0.5)
})

test_that("DirGLM accepts an explicit baseline prior", {
  dat <- data.frame(y = seq(0, 1, length.out = 6), x = -2:3)
  fit <- suppressWarnings(dirglm(
    y ~ x, dat, link = "logit",
    dirglmControl = dirglm.control(
      burnin = 0, thin = 1, save = 2, spt = sort(unique(dat$y)),
      dir_pr_parm = rep(1, 6), joint.update = FALSE, seed = 832
    )
  ))
  expect_true(all(is.finite(fit$samples$beta)))
  expect_true(fit$beta_acceptance >= 0 && fit$beta_acceptance <= 1)
})

test_that("an invalid exchangeable start uses the group-specific default", {
  dat <- update_test_data()
  dat$x <- dat$x - mean(dat$x)
  groups <- rep(1:2, each = 4)
  controls <- update_test_control()
  controls$corr <- "ex"
  controls$rhoStart <- -0.5
  controls$n_grid <- 71
  controls$r_grid <- 53
  expect_warning(bad_start <- cdpglm(
    y ~ x, dat, group_index = groups,
    cdpglmControl = do.call(cdpglm.control, controls)
  ), "rhoStart")
  controls$rhoStart <- NULL
  default_start <- cdpglm(
    y ~ x, dat, group_index = groups,
    cdpglmControl = do.call(cdpglm.control, controls)
  )
  expect_equal(bad_start$samples, default_start$samples)
})

test_that("rho proposals rounded to a boundary are rejected", {
  dat <- update_test_data()
  controls <- update_test_control()
  controls$corr <- "ar1"
  controls$rho_proposal_sd <- 1e6
  controls$n_grid <- 71
  controls$r_grid <- 53
  fit <- suppressWarnings(cdpglm(
    y ~ x, dat, group_index = dat$id,
    cdpglmControl = do.call(cdpglm.control, controls)
  ))
  expect_equal(fit$samples$rho, c(0, 0))
})

test_that("eps has no imposed upper bound in the CRM sampler", {
  expect_silent(crm_sampler(
    M = 2, u = c(1, 2), zstar = 0.5, nstar = 2, theta = c(-1, 1),
    alpha = 1, min_y = 0, max_y = 1, eps = 1e-3, h = 0,
    n_grid = 31, r_grid = 41
  ))
})

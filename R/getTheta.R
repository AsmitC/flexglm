## g function (logit transformation from appendix)
g <- function(mu, m, M) log(mu - m) - log(M - mu)

## Clamp values into the interior of (m, M) for stable transformed updates
clamp_interior <- function(x, m, M, eps) {
  pmin(pmax(x, m + eps), M - eps)
}

#' Control arguments for \eqn{\theta} update algorithm
#'
#' This function returns control arguments for the \eqn{\theta} update
#' algorithm.
#' Each argument has a default value, which will be used unless a different
#' value is provided by the user.
#'
#' @param eps Convergence threshold for theta updates. Convergence is evaluated
#' separately for each observation. An observation has converged when the
#' difference between \eqn{b'(\theta)} and \eqn{\mu} is less than
#' \code{epsTheta}.
#' @param maxiter Maximum number of iterations.
#' @param maxhalf Maximum number of half steps allowed per iteration if the
#' convergence criterion does not improve.
#' @param maxtheta Absolute value of theta is not allowed to exceed
#' \code{maxtheta}.
#' @param logit Logical for whether logit transformation should be used. Use of
#' this stabilizing transformation appears to be faster in general. Default is
#' TRUE.
#' @param logsumexp Logical argument for whether log-sum-exp trick should be
#' used. This may improve numerical stability at the expense of computational
#' time.
#'
#' @return Object of S3 class "thetaControl", which is a list of control
#' arguments.
#'
#' @export
theta.control <- function(eps = 1e-10, maxiter = 100, maxhalf = 20,
                          maxtheta = 500,
                          logit = TRUE, logsumexp = FALSE) {
  thetaControl <- as.list(environment())
  class(thetaControl) <- "thetaControl"
  thetaControl
}

#' getTheta
#' Updates theta. Vectorized but only updates observations that have not
#' converged.
#'
#' @param spt Support of the observed response variable. (This is the set of
#' unique values observed, not the set of all possible values.)
#' @param f0 Values of the baseline distribution corresponding to the values of
#' spt
#' @param mu The fitted mean for each observation. Note these values must lie
#' strictly within the range of the support.
#' @param sampprobs Matrix of sampling probabilities. The number of rows should
#' equal the number of observations, and the number of columns should equal the
#' number of unique observed support points.
#' @param ySptIndex Vector containing index of each obervation's response value
#' within the \code{spt} vector. This is only needed to calculate the
#' log-likelihood after each update.
#' @param thetaStart Vector of starting values. One value per observation. If
#' \code{NULL}, zero is used as the starting value for each observation.
#' @param thetaControl Object of class \code{thetaControl}, which is a list of
#' control arguments returned by the \code{thetaControl} function.
#'
#' @return List containing the following:
#' \itemize{
#' \item \code{theta} Updated values.
#' \item \code{fTilt} Matrix containing the exponentially tilted distribution
#' for each
#' observation, i.e. f(y|X=x). Each column corresponds to an observation and
#' sums to one.
#' \item \code{bPrime} Vector containing the mean of the exponentially tilted
#' distribution
#' for each observation. Should match \code{mu} argument very closely.
#' \item \code{bPrime2} Vector containing the variance of the exponentially
#' tilted
#' distribution for each observation.
#' \item \code{fTiltSW} Matrix containing the exponentially tilted distribution
#' for each
#' observation, conditional on that observation being sampled, i.e. f(y|X=x,
#' S=1).
#' If \code{sampprobs=NULL}, then \code{fTiltSW} matches \code{fTilt}.
#' \item \code{bPrimeSW} Vector containing the mean for each observation,
#' conditional
#' on that observation being sampled. If \code{sampprobs=NULL}, then
#' \code{bPrimeSW}
#' matches \code{bPrime}.
#' \item \code{bPrime2SW} Vector containing the variance for each observation,
#' conditional
#' on that observation being sampled. If \code{sampprobs=NULL}, then
#' \code{bPrime2SW}
#' matches \code{bPrime2}.
#' \item \code{llik} Semiparametric log-likelihood, evaluated at the current
#' beta
#' and f0 values. If sampling weights are used, then the log-likelihood is
#' conditional
#' on each observation being sampled.
#' \item \code{conv} Convergence indicator.
#' \item \code{iter} Number of iterations until convergence was reached.
#' }
#'
#' @keywords internal
getTheta <- function(spt, f0, mu, sampprobs = NULL, ySptIndex = NULL,
                     thetaStart = NULL, thetaControl = theta.control()) {
  if (!inherits(thetaControl, "thetaControl")) {
    stop(
      paste0("thetaControl must be an object of class 'thetaControl' ",
             "returned by theta.control()."),
      call. = FALSE
    )
  }

  getTheta_cpp(
    spt = spt,
    f0 = f0,
    mu = mu,
    sampprobs = sampprobs,
    ySptIndex = ySptIndex,
    thetaStart = thetaStart,
    eps = thetaControl$eps,
    maxiter = thetaControl$maxiter,
    maxhalf = thetaControl$maxhalf,
    maxtheta = thetaControl$maxtheta,
    logit = thetaControl$logit,
    logsumexp = thetaControl$logsumexp
  )
}

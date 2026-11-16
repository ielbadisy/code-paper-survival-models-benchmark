# R/dgp-rct-stress.R -- stressors layered on the trial DGP in R/dgp-rct.R.
# R/dgp-rct.R is shared with the letter and stays untouched; this file
# reuses rct_lp(), hazard_shapes and load_reference() and generalises the
# event-time model to
#   S(t | x, a) = exp(-(t / scale)^{k(x, a)} * exp(eta(x, a)))
# with
#   eta  = centred linear predictor (+ non-linear term if stress = "nle")
#   k    = arm-specific Weibull shape, multiplied by (1 + TVE_K * x_tve)
#          if stress = "tve" (a prognostic binary covariate whose hazard
#          ratio grows with time, beyond the arm-level shape difference)
# and censoring C ~ U[0, b * m(x)], with m(x) = 1 except for
# stress = "censdep" where m(x) = exp(-CENS_ALPHA * z(x)), z the
# standardised prognostic index: high-risk patients are censored earlier.
#
# stress in {"none", "nle", "tve", "censdep"}; "none" reproduces R/dgp-rct.R.

source(file.path(if (file.exists("R/dgp-rct.R")) "R" else "../R", "dgp-rct.R"))

TVE_K <- 0.5
CENS_ALPHA <- 0.7

# per dataset: covariate made non-linear (true effect on log scale), and the
# binary covariate given a time-varying effect
stress_spec <- list(
  pbc = list(nle_var = "bili", nle_beta = 0.85, tve_var = "hepato"),
  prostate = list(nle_var = "ap", nle_beta = 0.35, tve_var = "hx")
)

.stress_eta <- function(X, trt, ref, beta_trt, stress) {
  lp <- rct_lp(X, trt, ref, beta_trt)
  if (stress == "nle") {
    s <- stress_spec[[ref$name]]; v <- s$nle_var
    mu_lin <- mean(ref$X[[v]]); mu_log <- mean(log(ref$X[[v]]))
    # swap the linear term for a log term with the same reference centring
    lp <- lp - ref$beta[[v]] * (X[[v]] - mu_lin) + s$nle_beta * (log(X[[v]]) - mu_log)
  }
  lp
}

.stress_shape <- function(X, trt, ref, hazard, stress) {
  k <- hazard_shapes[[hazard]][trt + 1]
  if (stress == "tve") k <- k * (1 + TVE_K * X[[stress_spec[[ref$name]]$tve_var]])
  k
}

stress_true_surv <- function(X, trt, ref, beta_trt, hazard, stress, times) {
  eta <- exp(.stress_eta(X, trt, ref, beta_trt, stress))
  k <- .stress_shape(X, trt, ref, hazard, stress)
  outer(seq_along(eta), times, function(i, t) exp(-(t / ref$scale)^k[i] * eta[i]))
}

stress_draw_times <- function(X, trt, ref, beta_trt, hazard, stress) {
  eta <- exp(.stress_eta(X, trt, ref, beta_trt, stress))
  k <- .stress_shape(X, trt, ref, hazard, stress)
  ref$scale * (-log(runif(length(eta))) / eta)^(1 / k)
}

.cens_mult <- function(X, ref, stress) {
  if (stress != "censdep") return(rep(1, nrow(X)))
  pi0 <- rct_lp(ref$X, rep(0L, nrow(ref$X)), ref, 0)   # reference prognostic index
  z <- (rct_lp(X, rep(0L, nrow(X)), ref, 0) - mean(pi0)) / stats::sd(pi0)
  exp(-CENS_ALPHA * z)
}

calibrate_stress_censoring <- function(ref, beta_trt, hazard, stress, target,
                                       n_pilot = 1e5, seed = 1) {
  set.seed(seed)
  X <- ref$X[sample.int(nrow(ref$X), n_pilot, TRUE), , drop = FALSE]
  trt <- rbinom(n_pilot, 1, 0.5)
  T <- stress_draw_times(X, trt, ref, beta_trt, hazard, stress)
  m <- .cens_mult(X, ref, stress)
  # P(C < T) with C ~ U[0, b m] is mean(pmin(T, b m) / (b m))
  f <- function(b) mean(pmin(T, b * m) / (b * m)) - target
  uniroot(f, c(1e-6, 1e3) * stats::quantile(T, 0.99), tol = 1e-8)$root
}

simulate_stress <- function(ref, n, beta_trt, hazard, stress, cens_b) {
  X <- ref$X[sample.int(nrow(ref$X), n, TRUE), , drop = FALSE]; rownames(X) <- NULL
  trt <- rbinom(n, 1, 0.5)
  T <- stress_draw_times(X, trt, ref, beta_trt, hazard, stress)
  C <- runif(n, 0, cens_b * .cens_mult(X, ref, stress))
  d <- data.frame(time = pmin(T, C), status = as.integer(T <= C), trt = trt, X)
  attr(d, "truth") <- list(X = X, trt = trt)
  d
}

# R/dgp.R -- trial-shaped DGPs on the two reference RCTs of Graf, Todd & Baksh (2026, BIMJ 68:e70171):
# PBC (no treatment-covariate interactions) and VACURG prostate (three interactions).
#   S(t | x, a) = exp(-(t / scale)^k(x, a) * exp(eta(x, a))), C ~ U[0, b * m(x)]
# Covariates are resampled trial rows (plasmode), treatment re-randomised 1:1; eta uses covariates
# centred at the trial means, so beta_trt is the log-HR at the average patient. Stressors:
#   nle: one covariate enters as log(x); tve: shape k multiplied by 1 + TVE_K * x_tve (binary x);
#   censdep: m(x) = exp(-CENS_ALPHA * z(x)), z the standardised prognostic index.
suppressMessages(library(survival))
TVE_K <- 0.5; CENS_ALPHA <- 0.7
DATA_REF <- if (dir.exists("data-ref")) "data-ref" else "../data-ref"   # repo root or an arm folder
# Weibull shapes (control, treated); "nph" gives crossing hazards
HAZ <- list(dec = c(0.8, 0.8), const = c(1, 1), inc = c(2, 2), nph = c(2, 5))
STRESS <- list(pbc = list(nle_var = "bili", nle_beta = 0.85, tve_var = "hepato"),
               prostate = list(nle_var = "ap", nle_beta = 0.35, tve_var = "hx"))

load_reference <- function(dataset) {
  if (dataset == "pbc") {
    d <- subset(survival::pbc, !is.na(trt))
    for (v in c("chol", "copper", "trig", "platelet")) d[[v]][is.na(d[[v]])] <- mean(d[[v]], na.rm = TRUE)
    d$sex <- as.integer(d$sex == "f")
    covs <- c("age", "sex", "ascites", "hepato", "spiders", "edema", "bili", "chol", "albumin", "copper",
              "alk.phos", "ast", "trig", "platelet", "protime", "stage")
    # Graf et al. Section 2.2.2, Z2..Z17
    beta <- c(0.026, -0.218, 0.338, 0.227, 0.071, 0.481, 0.086, 0.0004, -0.799, 0.003, -0.00002, 0.004,
              -0.002, 0.0002, 0.276, 0.365)
    list(name = "pbc", X = d[covs], beta = setNames(beta, covs), inter = NULL, scale = 2241.74)
  } else {
    d <- read.csv(file.path(DATA_REF, "adv_prostate_ca/adv_prostate_ca.csv"))
    covs <- c("age", "wt", "sbp", "dbp", "sz", "ap", "hg", "sg", "pf", "hx", "bm", "ekg")
    beta <- c(-0.006, -0.01, -0.016, 0.02, 0.014, 0.0001, -0.006, 0.074, 0.333, 0.467, 0.63, 0.316)
    list(name = "prostate", X = d[covs], beta = setNames(beta, covs),
         inter = c(age = 0.059, bm = -0.612, ap = -0.0003), scale = 39.2)
  }
}

rct_lp <- function(X, trt, ref, beta_trt) {
  Xc <- sweep(as.matrix(X), 2, colMeans(ref$X))
  lp <- drop(Xc %*% ref$beta[colnames(Xc)]) + beta_trt * trt
  for (v in names(ref$inter)) lp <- lp + ref$inter[[v]] * trt * Xc[, v]
  lp
}
.stress_eta <- function(X, trt, ref, beta_trt, stress) {
  lp <- rct_lp(X, trt, ref, beta_trt)
  if (stress != "nle") return(lp)
  s <- STRESS[[ref$name]]; v <- s$nle_var
  # swap the linear term for a log term with the same reference centring
  lp - ref$beta[[v]] * (X[[v]] - mean(ref$X[[v]])) + s$nle_beta * (log(X[[v]]) - mean(log(ref$X[[v]])))
}
.stress_shape <- function(X, trt, ref, hazard, stress) {
  k <- HAZ[[hazard]][trt + 1]
  if (stress == "tve") k * (1 + TVE_K * X[[STRESS[[ref$name]]$tve_var]]) else k
}
stress_true_surv <- function(X, trt, ref, beta_trt, hazard, stress, times) {
  eta <- exp(.stress_eta(X, trt, ref, beta_trt, stress)); k <- .stress_shape(X, trt, ref, hazard, stress)
  outer(seq_along(eta), times, function(i, t) exp(-(t / ref$scale)^k[i] * eta[i]))
}
stress_draw_times <- function(X, trt, ref, beta_trt, hazard, stress) {
  eta <- exp(.stress_eta(X, trt, ref, beta_trt, stress)); k <- .stress_shape(X, trt, ref, hazard, stress)
  ref$scale * (-log(runif(length(eta))) / eta)^(1 / k)
}
.cens_mult <- function(X, ref, stress) {
  if (stress != "censdep") return(rep(1, nrow(X)))
  pi0 <- rct_lp(ref$X, rep(0L, nrow(ref$X)), ref, 0)
  exp(-CENS_ALPHA * (rct_lp(X, rep(0L, nrow(X)), ref, 0) - mean(pi0)) / sd(pi0))
}
# upper bound b giving the target censored share: P(C < T) = mean(pmin(T, b m) / (b m))
calibrate_stress_censoring <- function(ref, beta_trt, hazard, stress, target, n_pilot = 1e5, seed = 1) {
  set.seed(seed)
  X <- ref$X[sample.int(nrow(ref$X), n_pilot, TRUE), , drop = FALSE]
  T <- stress_draw_times(X, rbinom(n_pilot, 1, 0.5), ref, beta_trt, hazard, stress)
  m <- .cens_mult(X, ref, stress)
  uniroot(function(b) mean(pmin(T, b * m) / (b * m)) - target, c(1e-6, 1e3) * quantile(T, 0.99), tol = 1e-8)$root
}
simulate_stress <- function(ref, n, beta_trt, hazard, stress, cens_b) {
  X <- ref$X[sample.int(nrow(ref$X), n, TRUE), , drop = FALSE]; rownames(X) <- NULL
  trt <- rbinom(n, 1, 0.5)
  T <- stress_draw_times(X, trt, ref, beta_trt, hazard, stress)
  C <- runif(n, 0, cens_b * .cens_mult(X, ref, stress))
  structure(data.frame(time = pmin(T, C), status = as.integer(T <= C), trt = trt, X),
            truth = list(X = X, trt = trt))
}

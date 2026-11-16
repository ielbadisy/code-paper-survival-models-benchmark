# R/dgp-rct.R -- trial-shaped DGPs built on the two reference RCTs used by
# Graf, Todd & Baksh (2026, Biometrical Journal 68:e70171): PBC (no
# treatment-covariate interactions) and VACURG advanced prostate cancer
# (three treatment-covariate interactions). Shared by letter/ and the paper.
#
# Differences from Graf et al. (documented, deliberate):
# * covariates are resampled rows of the real trial (plasmode), not a fitted
#   copula; treatment is re-randomised 1:1. This keeps the joint covariate
#   distribution exactly.
# * the linear predictor uses covariates centred at the reference-data means,
#   so beta_trt is the treatment log-HR at the average patient. With the
#   published (uncentred) prostate coefficients, rx:age = 0.059 makes the
#   effective treatment log-HR at age 72 about beta_trt + 4.2.
#   centre = FALSE reproduces the uncentred version.
# * every simulated dataset carries the TRUE survival function, so
#   calibration can be measured against the truth.

suppressMessages(library(survival))

load_reference <- function(dataset = c("pbc", "prostate")) {
  dataset <- match.arg(dataset)
  if (dataset == "pbc") {
    d <- subset(survival::pbc, !is.na(trt))
    for (v in c("chol", "copper", "trig", "platelet")) d[[v]][is.na(d[[v]])] <- mean(d[[v]], na.rm = TRUE)
    d$sex <- as.integer(d$sex == "f")
    d$trt <- as.integer(d$trt == 1)  # 1 = D-penicillamine
    covs <- c("age", "sex", "ascites", "hepato", "spiders", "edema", "bili", "chol",
              "albumin", "copper", "alk.phos", "ast", "trig", "platelet", "protime", "stage")
    # Graf et al. Section 2.2.2, Z2..Z17, in survival::pbc column order
    beta <- c(0.026, -0.218, 0.338, 0.227, 0.071, 0.481, 0.086, 0.0004, -0.799, 0.003,
              -0.00002, 0.004, -0.002, 0.0002, 0.276, 0.365)
    list(name = "pbc", X = d[covs], beta = setNames(beta, covs), inter = NULL,
         scale = 2241.74, trt = "trt")
  } else {
    f <- file.path("data-ref", "adv_prostate_ca", "adv_prostate_ca.csv")
    if (!file.exists(f)) f <- file.path("..", f)
    d <- read.csv(f)
    covs <- c("age", "wt", "sbp", "dbp", "sz", "ap", "hg", "sg", "pf", "hx", "bm", "ekg")
    beta <- c(-0.006, -0.01, -0.016, 0.02, 0.014, 0.0001, -0.006, 0.074, 0.333, 0.467, 0.63, 0.316)
    list(name = "prostate", X = d[covs], beta = setNames(beta, covs),
         inter = c(age = 0.059, bm = -0.612, ap = -0.0003), scale = 39.2, trt = "rx")
  }
}

# hazard: "dec" (gamma 0.8), "const" (1), "inc" (2, as in Graf's text; their
# figures say 6), "inc6" (6), "nph" (gamma0 = 2 control, gamma1 = 5 treated)
hazard_shapes <- list(dec = c(0.8, 0.8), const = c(1, 1), inc = c(2, 2),
                      inc6 = c(6, 6), nph = c(2, 5))

rct_lp <- function(X, trt, ref, beta_trt, centre = TRUE) {
  mu <- if (centre) colMeans(ref$X) else setNames(rep(0, ncol(ref$X)), names(ref$X))
  Xc <- sweep(as.matrix(X), 2, mu)
  lp <- drop(Xc %*% ref$beta[colnames(Xc)]) + beta_trt * trt
  if (!is.null(ref$inter))
    for (v in names(ref$inter)) lp <- lp + ref$inter[[v]] * trt * Xc[, v]
  lp
}

# true S(t | x, arm) for each row of X (rows) at each time (cols)
rct_true_surv <- function(X, trt, ref, beta_trt, hazard, times, centre = TRUE) {
  g <- hazard_shapes[[hazard]][trt + 1]
  eta <- exp(rct_lp(X, trt, ref, beta_trt, centre))
  outer(seq_along(eta), times, function(i, t) exp(-(t / ref$scale)^g[i] * eta[i]))
}

draw_rct_times <- function(X, trt, ref, beta_trt, hazard, centre = TRUE) {
  g <- hazard_shapes[[hazard]][trt + 1]
  eta <- exp(rct_lp(X, trt, ref, beta_trt, centre))
  ref$scale * (-log(runif(length(eta))) / eta)^(1 / g)
}

# upper bound b of U[0, b] censoring giving the target censoring proportion
calibrate_rct_censoring <- function(ref, beta_trt, hazard, target, centre = TRUE,
                                    n_pilot = 1e5, seed = 1) {
  set.seed(seed)
  idx <- sample.int(nrow(ref$X), n_pilot, TRUE)
  trt <- rbinom(n_pilot, 1, 0.5)
  T <- draw_rct_times(ref$X[idx, ], trt, ref, beta_trt, hazard, centre)
  # P(C < T) with C ~ U[0, b] is mean(pmin(T, b) / b)
  f <- function(b) mean(pmin(T, b) / b) - target
  uniroot(f, c(1e-6, 1e3) * stats::quantile(T, 0.99), tol = 1e-8)$root
}

simulate_rct <- function(ref, n, beta_trt, hazard, cens_b, centre = TRUE) {
  idx <- sample.int(nrow(ref$X), n, TRUE)
  X <- ref$X[idx, , drop = FALSE]; rownames(X) <- NULL
  trt <- rbinom(n, 1, 0.5)
  T <- draw_rct_times(X, trt, ref, beta_trt, hazard, centre)
  C <- runif(n, 0, cens_b)
  d <- data.frame(time = pmin(T, C), status = as.integer(T <= C), trt = trt, X)
  attr(d, "truth") <- list(X = X, trt = trt)
  d
}

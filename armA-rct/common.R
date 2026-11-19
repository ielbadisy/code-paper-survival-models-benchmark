# armA-rct/common.R -- design, learners and scoring shared by tune-run.R and
# sim-run.R (v2, 2026-11-17). Sourced from armA-rct/.

suppressMessages({ library(survalis); library(survival); library(mcstatsim) })
source("../R/dgp-rct-stress.R"); source("../R/benefit.R"); source("../R/sim-core.R")
options(rf.cores = 1L, mc.cores = 1L)

N_TEST <- 500L
REFS <- list(pbc = load_reference("pbc"), prostate = load_reference("prostate"))
LEARNERS <- c("coxph", "cox_spec", "glmnet", "coxboost", "rpsurv", "rpsurv_tve", "nbsurv", "fastgbm",
              "densemlp", "bart", "xgboost", "rsf", "orsf", "blackboost")
# fixed (untuned) arguments: single thread, capped BART chain; densemlp
# epochs only caps early stopping
BASE <- list(bart = list(ntree = 50, ndpost = 200, nskip = 100, printevery = 1e6),
             fastgbm = list(threads = 1L), densemlp = list(ncores = 1L, epochs = 200, weight_decay = 1e-3),
             orsf = list(n_tree = 300, n_thread = 1L), rsf = list(ntree = 300))
DW_RANGE <- "range:0.2-0.8"   # decision-weighted loss over risk thresholds at t*
PKGS <- c("survalis", "nbsurv", "fastgbm", "densemlp", "mcstatsim", "BART", "xgboost",
          "randomForestSRC", "aorsf", "mboost", "CoxBoost", "glmnet", "rpsurv")

base_grid <- expand.grid(dataset = c("pbc", "prostate"), hazard = c("dec", "const", "inc", "nph"),
                         beta_trt = c(0, -0.4, 0.8), n = c(100L, 200L, 400L), cens = c(0.3, 0.6),
                         stress = "none", stringsAsFactors = FALSE)
stress_grid <- expand.grid(dataset = c("pbc", "prostate"), hazard = c("const", "nph"),
                           beta_trt = -0.4, n = 200L, cens = 0.3,
                           stress = c("nle", "tve", "censdep"), stringsAsFactors = FALSE)
GRID <- rbind(base_grid, stress_grid)

cell_key <- function(ds, hz, bt, st, cs) paste(ds, hz, bt, st, cs)
# deciles 10..90 of the true event times restricted to the follow-up window
# [0, tau], tau = 90th percentile of min(T, C): keeps every time inside the
# censoring support so IPCW metrics exist; column 5 is t*
eval_times <- function(ref, beta_trt, hazard, stress, cens_b) {
  set.seed(7)
  X <- ref$X[sample.int(nrow(ref$X), 1e5, TRUE), , drop = FALSE]
  T <- stress_draw_times(X, rbinom(1e5, 1, 0.5), ref, beta_trt, hazard, stress)
  C <- runif(1e5, 0, cens_b * .cens_mult(X, ref, stress))
  tau <- unname(quantile(pmin(T, C), 0.9))
  unname(quantile(T[T <= tau], seq(0.1, 0.9, 0.1)))
}
prepare_grid <- function(g) {
  g$cens_b <- mapply(function(ds, hz, bt, st, cs) calibrate_stress_censoring(REFS[[ds]], bt, hz, st, cs, n_pilot = 5e4),
                     g$dataset, g$hazard, g$beta_trt, g$stress, g$cens)
  g
}
grid_times <- function(g) {
  k <- unique(g[c("dataset", "hazard", "beta_trt", "stress", "cens", "cens_b")])
  setNames(lapply(seq_len(nrow(k)), function(i) with(k[i, ], eval_times(REFS[[dataset]], beta_trt, hazard, stress, cens_b))),
           with(k, cell_key(dataset, hazard, beta_trt, stress, cens)))
}

# oracle-form Cox: strata(trt) + the true treatment-covariate interactions
fit_cox_spec <- function(train, ref) {
  ix <- names(ref$inter)
  for (v in ix) train[[paste0("trtx_", v)]] <- train$trt * train[[v]]
  rhs <- paste(c(names(ref$X), if (length(ix)) paste0("trtx_", ix)), collapse = " + ")
  structure(list(fit = coxph(as.formula(paste("Surv(time, status) ~", rhs, "+ strata(trt)")), train), ix = ix),
            class = "coxspec")
}
# S(t | x) = exp(-H0_arm(t) exp(x'b)), Breslow H0 per stratum (survfit with
# strata on 500 rows takes about a minute)
predict_cox_spec <- function(obj, newdata, times) {
  for (v in obj$ix) newdata[[paste0("trtx_", v)]] <- newdata$trt * newdata[[v]]
  b <- coef(obj$fit); b[is.na(b)] <- 0
  eta <- exp(drop(as.matrix(newdata[names(b)]) %*% b))
  bh <- basehaz(obj$fit, centered = FALSE)
  S <- matrix(NA_real_, nrow(newdata), length(times))
  for (a in 0:1) {
    i <- which(newdata$trt == a); if (!length(i)) next
    h <- bh[bh$strata == paste0("trt=", a), ]
    S[i, ] <- exp(-outer(eta[i], stats::stepfun(h$time, c(0, h$hazard))(times)))
  }
  S
}

# CoxBoost with trt unpenalised (survalis does not expose unpen.index);
# penalty = pmult x events (CoxBoost default pmult 9)
covariate_matrix <- function(d) as.matrix(d[setdiff(names(d), c("time", "status"))])
fit_coxboost_unpen <- function(train, hp) {
  x <- covariate_matrix(train)
  a <- modifyList(list(stepno = 100, pmult = 9), hp)
  m <- CoxBoost::CoxBoost(time = train$time, status = train$status, x = x,
                          unpen.index = which(colnames(x) == "trt"), stepno = a$stepno,
                          penalty = a$pmult * sum(train$status == 1))
  structure(list(model = m, cols = colnames(x)), class = "coxboostunpen")
}

# trt is unpenalised in glmnet and coxboost so shrinkage cannot remove the
# treatment effect; stdout of chatty engines (BART) is swallowed
fit_any <- function(learner, train, ref, hp = NULL) {
  fo <- Surv(time, status) ~ .
  hp <- modifyList(if (is.null(BASE[[learner]])) list() else BASE[[learner]], if (is.null(hp)) list() else hp)
  pf <- as.numeric(setdiff(names(train), c("time", "status")) != "trt")
  utils::capture.output(f <- switch(learner,
    cox_spec = fit_cox_spec(train, ref),
    rpsurv_tve = fit(fo, train, model = "rpsurv", spec = list(tve = "trt")),
    coxboost = fit_coxboost_unpen(train, hp),
    glmnet = fit(fo, train, model = "glmnet", spec = c(hp, list(penalty.factor = pf))),
    fit(fo, train, model = learner, spec = hp)))
  f
}
predict_any <- function(obj, newdata, times) {
  if (inherits(obj, "coxspec")) return(predict_cox_spec(obj, newdata, times))
  if (inherits(obj, "coxboostunpen"))
    return(matrix(predict(obj$model, newdata = covariate_matrix(newdata)[, obj$cols, drop = FALSE],
                          times = times, type = "risk"), nrow(newdata), length(times)))
  as.matrix(predict(obj, newdata = newdata, times = times))
}
pkg_versions <- function() vapply(PKGS, function(p) as.character(packageVersion(p)), "")

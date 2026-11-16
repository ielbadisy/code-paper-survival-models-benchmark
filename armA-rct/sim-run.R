# armA-rct/sim-run.R -- Arm A: trial-shaped benchmark built on the two
# reference RCTs of Graf, Todd & Baksh (2026, BIMJ 68:e70171).
# Every replicate: simulate a training trial and an independent test trial
# (N = 500), fit each learner on (trt + covariates), predict S(t | x) on the
# test set with the observed arm and with trt set to 1 and to 0, then score
#   * prediction: Harrell C, Uno C, td-AUC, IBS, Brier/ECE/survacc at t*
#   * truth calibration: mean |S_hat - S_true| over the grid and at t*
#   * benefit at t*: RMSE, bias, calibration slope, c-for-benefit,
#     policy regret and gain (R/benefit.R)
# t* = median of the true event-time distribution of the cell, restricted to
# the follow-up window (see eval_times).
# Learner settings: tuned once per dataset x hazard by tune-run.R
# (results/tuned_armA.rds); trt is unpenalised in glmnet and coxboost so
# that shrinkage cannot remove the treatment effect (PROTOCOL.md Section 10).
# Results are checkpointed by runsim (save_path) every CKPT_EVERY replicates;
# rerunning the same command resumes from the checkpoint.
# Run from armA-rct/:  Rscript sim-run.R <reps> <cores> [smoke|pilot]

suppressMessages({ library(survalis); library(survival); library(mcstatsim) })
source("../R/dgp-rct-stress.R")
source("../R/benefit.R")
source("../R/sim-core.R")

N_TEST <- 500L
REFS <- list(pbc = load_reference("pbc"), prostate = load_reference("prostate"))

# learners: survalis set minus exclusions (see PROTOCOL.md, Learners), plus
# two hand-built comparators: cox_spec (oracle-form Cox: strata(trt) + the
# true treatment-covariate interactions) and rpsurv_tve (Royston-Parmar
# with a time-varying treatment effect)
# Held out pending capped specs (default fits on PBC n = 400, 2 threads):
# xgboost 172 s (C 0.645, IBS 0.203), densemlp 385 s, bart 81 s.
LEARNERS <- c("coxph", "glmnet", "rpsurv", "rpsurv_tve", "cox_spec", "flexsurvreg",
              "coxboost", "fastgbm", "rsf", "ranger", "orsf", "cforest",
              "blackboost", "nbsurv", "survsvm")
HELD_OUT <- c("xgboost", "densemlp", "bart")   # future work (cost), PROTOCOL.md Section 10
TUNED <- if (file.exists("results/tuned_armA.rds")) readRDS("results/tuned_armA.rds") else list()
CKPT_EVERY <- 5L

# ---- grid -----------------------------------------------------------------
base_grid <- expand.grid(dataset = c("pbc", "prostate"), hazard = c("dec", "const", "inc", "nph"),
                         beta_trt = c(0, -0.4, 0.8), n = c(100L, 200L, 400L), cens = c(0.3, 0.6),
                         stress = "none", stringsAsFactors = FALSE)
stress_grid <- expand.grid(dataset = c("pbc", "prostate"), hazard = c("const", "nph"),
                           beta_trt = -0.4, n = 200L, cens = 0.3,
                           stress = c("nle", "tve", "censdep"), stringsAsFactors = FALSE)
grid <- rbind(base_grid, stress_grid)

# ---- per-cell constants: censoring bound and evaluation times ------------
# key includes the censoring level: the grid depends on the follow-up window
cell_key <- function(ds, hz, bt, st, cs) paste(ds, hz, bt, st, cs)
# evaluation grid: deciles 10..90 of the true event-time distribution
# restricted to the follow-up window [0, tau], tau = 90th percentile of the
# observed follow-up min(T, C). Without the restriction the upper times fall
# beyond the censoring support U[0, b] (IPCW undefined: IBS/td-AUC all NA in
# the 2026-11-15 pilot, and t* > b in some 60% censoring cells).
eval_times <- function(ref, beta_trt, hazard, stress, cens_b) {
  set.seed(7)
  X <- ref$X[sample.int(nrow(ref$X), 1e5, TRUE), , drop = FALSE]
  T <- stress_draw_times(X, rbinom(1e5, 1, 0.5), ref, beta_trt, hazard, stress)
  C <- runif(1e5, 0, cens_b * .cens_mult(X, ref, stress))
  tau <- unname(quantile(pmin(T, C), 0.9))
  unname(quantile(T[T <= tau], seq(0.1, 0.9, 0.1)))   # 9 times; column 5 is t*
}
prepare_grid <- function(grid) {
  grid$cens_b <- mapply(function(ds, hz, bt, st, cs)
    calibrate_stress_censoring(REFS[[ds]], bt, hz, st, cs, n_pilot = 5e4),
    grid$dataset, grid$hazard, grid$beta_trt, grid$stress, grid$cens)
  grid
}

# ---- hand-built comparators ----------------------------------------------
fit_cox_spec <- function(train, ref) {
  ix <- names(ref$inter)
  for (v in ix) train[[paste0("trtx_", v)]] <- train$trt * train[[v]]
  rhs <- paste(c(names(ref$X), if (length(ix)) paste0("trtx_", ix)), collapse = " + ")
  structure(list(fit = coxph(as.formula(paste("Surv(time, status) ~", rhs, "+ strata(trt)")), train),
                 ix = ix), class = "coxspec")
}
# S(t | x) = exp(-H0_arm(t) exp(x'b)), Breslow H0 per stratum; vectorised
# (survfit.coxph with 500 newdata rows and strata takes about a minute)
predict_cox_spec <- function(obj, newdata, times) {
  for (v in obj$ix) newdata[[paste0("trtx_", v)]] <- newdata$trt * newdata[[v]]
  b <- coef(obj$fit); b[is.na(b)] <- 0
  eta <- exp(drop(as.matrix(newdata[names(b)]) %*% b))
  bh <- basehaz(obj$fit, centered = FALSE)
  S <- matrix(NA_real_, nrow(newdata), length(times))
  for (a in 0:1) {
    i <- which(newdata$trt == a)
    if (!length(i)) next
    h <- bh[bh$strata == paste0("trt=", a), ]
    H0 <- stats::stepfun(h$time, c(0, h$hazard))(times)
    S[i, ] <- exp(-outer(eta[i], H0))
  }
  S
}

# CoxBoost with trt unpenalised (unpen.index); survalis fit() does not
# expose unpen.index. Defaults as survalis: stepno 100, penalty 9 x events.
covariate_matrix <- function(d) as.matrix(d[setdiff(names(d), c("time", "status"))])
fit_coxboost_unpen <- function(train, hp = NULL) {
  x <- covariate_matrix(train)
  a <- list(stepno = 100, penalty = 9 * sum(train$status == 1), criterion = "pscore",
            stepsize.factor = 1)
  a[names(hp)] <- hp
  m <- CoxBoost::CoxBoost(time = train$time, status = train$status, x = x,
                          unpen.index = which(colnames(x) == "trt"),
                          stepno = a$stepno, penalty = a$penalty, criterion = a$criterion,
                          stepsize.factor = a$stepsize.factor)
  structure(list(model = m, cols = colnames(x)), class = "coxboostunpen")
}

fit_any <- function(learner, train, ref, hp = NULL) {
  fo <- Surv(time, status) ~ .
  if (is.null(hp)) hp <- list()
  switch(learner,
    cox_spec = fit_cox_spec(train, ref),
    rpsurv_tve = fit(fo, train, model = "rpsurv", spec = list(tve = "trt")),
    coxboost = fit_coxboost_unpen(train, hp),
    glmnet = fit(fo, train, model = "glmnet", spec = c(hp, list(penalty.factor =
      as.numeric(setdiff(names(train), c("time", "status")) != "trt")))),
    fit(fo, train, model = learner, spec = hp))
}
predict_any <- function(obj, newdata, times) {
  if (inherits(obj, "coxspec")) predict_cox_spec(obj, newdata, times)
  else if (inherits(obj, "coxboostunpen"))
    matrix(predict(obj$model, newdata = covariate_matrix(newdata)[, obj$cols, drop = FALSE],
                   times = times, type = "risk"), nrow(newdata), length(times))
  else as.matrix(predict(obj, newdata = newdata, times = times))
}

# ---- one replicate --------------------------------------------------------
pred_names <- c(metric_names, truth_metric_names)
ben_names <- c("benefit_rmse", "benefit_bias", "benefit_cal_int", "benefit_cal_slope",
               "benefit_cindex", "policy_regret", "policy_gain")
all_names <- c(pred_names, ben_names, "secs")

sim_one <- function(dataset, hazard, beta_trt, n, cens, stress, cens_b) {
  ref <- REFS[[dataset]]
  times <- TIMES[[cell_key(dataset, hazard, beta_trt, stress, cens)]]
  train <- simulate_stress(ref, n, beta_trt, hazard, stress, cens_b)
  test <- simulate_stress(ref, N_TEST, beta_trt, hazard, stress, cens_b)
  attr(train, "truth") <- NULL
  tr <- attr(test, "truth"); attr(test, "truth") <- NULL
  S_true <- stress_true_surv(tr$X, tr$trt, ref, beta_trt, hazard, stress, times)
  S1 <- stress_true_surv(tr$X, rep(1L, N_TEST), ref, beta_trt, hazard, stress, times[5])[, 1]
  S0 <- stress_true_surv(tr$X, rep(0L, N_TEST), ref, beta_trt, hazard, stress, times[5])[, 1]
  y <- Surv(test$time, test$status)
  d1 <- test; d1$trt <- 1L; d0 <- test; d0$trt <- 0L

  do.call(rbind, lapply(LEARNERS, function(l) {
    t0 <- proc.time()[["elapsed"]]
    v <- tryCatch({
      f <- fit_any(l, train, ref, TUNED[[dataset]][[hazard]][[l]])
      S <- predict_any(f, test, times)
      bhat <- predict_any(f, d1, times[5])[, 1] - predict_any(f, d0, times[5])[, 1]
      if (anyNA(S) || anyNA(bhat)) stop("NA predictions")
      c(score_matrix(S, y, times, 5), score_truth(S, S_true, y, times, 5),
        benefitmetrics(bhat, S1, S0))
    }, error = function(e) setNames(rep(NA_real_, length(all_names) - 1), all_names[-length(all_names)]))
    v["secs"] <- proc.time()[["elapsed"]] - t0
    data.frame(learner = l, metric = names(v), value = unname(v))
  }))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  reps <- if (length(args) >= 1) as.integer(args[1]) else 200L
  cores <- if (length(args) >= 2) as.integer(args[2]) else 2L
  smoke <- length(args) >= 3 && args[3] == "smoke"
  pilot <- length(args) >= 3 && args[3] == "pilot"
  Sys.setenv(OMP_NUM_THREADS = 1)
  if (smoke) grid <- grid[c(which(grid$dataset == "pbc" & grid$hazard == "nph" & grid$beta_trt == -0.4 &
                                  grid$n == 200 & grid$cens == 0.3 & grid$stress == "none"),
                            which(grid$dataset == "prostate" & grid$hazard == "nph" & grid$beta_trt == -0.4 &
                                  grid$n == 200 & grid$cens == 0.3 & grid$stress == "none"),
                            which(grid$dataset == "prostate" & grid$stress == "tve" & grid$hazard == "nph")), ]
  # pilot: every dataset x hazard at the two extremes (n 100, 60% censoring,
  # beta 0.8; n 400, 30% censoring, beta -0.4) plus all 12 stress cells
  if (pilot) grid <- grid[(grid$stress == "none" & ((grid$n == 100 & grid$cens == 0.6 & grid$beta_trt == 0.8) |
                                                    (grid$n == 400 & grid$cens == 0.3 & grid$beta_trt == -0.4))) |
                          grid$stress != "none", ]
  rownames(grid) <- NULL
  if (!length(TUNED)) message("results/tuned_armA.rds not found: all learners at defaults")
  grid <- prepare_grid(grid)
  TIMES <- list()
  for (i in seq_len(nrow(grid))) with(grid[i, ], {
    k <- cell_key(dataset, hazard, beta_trt, stress, cens)
    if (is.null(TIMES[[k]])) TIMES[[k]] <<- eval_times(REFS[[dataset]], beta_trt, hazard, stress, cens_b)
  })
  dir.create("results", showWarnings = FALSE)
  tag <- if (smoke) "smoke" else if (pilot) sprintf("pilot%d", reps) else sprintf("%dreps", reps)
  saveRDS(list(grid = grid, times = TIMES, learners = LEARNERS, tuned = TUNED),
          sprintf("results/design_%s.rds", tag))
  raw <- runsim(n = reps, grid_params = grid, sim_func = sim_one, seed = 2027,
                num_cores = cores, backend = "fork", on_error = "warn",
                save_path = sprintf("results/checkpoint_%s.rds", tag), checkpoint_every = CKPT_EVERY)
  saveRDS(raw, sprintf("results/sim_raw_%s.rds", tag))
  write.csv(raw, sprintf("results/sim_raw_%s.csv", tag), row.names = FALSE)
}

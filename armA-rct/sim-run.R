# armA-rct/sim-run.R -- Arm A v2 evaluation (replicates 11-100), trial-shaped
# benchmark on the two reference RCTs of Graf, Todd & Baksh (2026, BIMJ).
# Every replicate: simulate a training trial (n) and an independent test trial
# (N = 500), fit each of the 14 learners with the settings frozen by
# tune-run.R (results/tuned_<ttag>.rds, keyed dataset x hazard x n), predict
# S(t | x) on the test set with the observed arm and with trt set to 1 and 0:
#   * prediction: Harrell C, td-AUC, IBS, IAE, ISE; Brier, ECE, survacc and
#     decision-weighted loss over risk 0.2-0.8 at t*
#   * truth: Uno C, mean |S_hat - S_true| (grid and t*), calibration slope of
#     S_true on S_hat at t*
#   * benefit at t*: RMSE, bias, calibration int/slope, c-for-benefit,
#     policy regret and gain (R/benefit.R)
# Inputs: common.R, results/tuned_<ttag>.rds
# Outputs: results/{design,checkpoint,sim_raw}_<tag>.rds (+ sim_raw csv)
# Run from armA-rct/:
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript sim-run.R <reps> <cores> [ttag]

source("common.R")

PRED <- c(metric_names, "dwrange", truth_metric_names, "cal_slope_true")
BEN <- c("benefit_rmse", "benefit_bias", "benefit_cal_int", "benefit_cal_slope", "benefit_cindex",
         "policy_regret", "policy_gain")
NA_ROW <- setNames(rep(NA_real_, length(PRED) + length(BEN)), c(PRED, BEN))

score_all <- function(S, S_true, y, times, bhat, S1, S0) {
  dw <- evaluate(S[, 5, drop = FALSE], y = y, times = times[5], metrics = DW_RANGE)$summary$mean
  cs <- coef(lm(S_true[, 5] ~ S[, 5]))[2]
  c(score_matrix(S, y, times, 5), dwrange = dw, score_truth(S, S_true, y, times, 5),
    cal_slope_true = unname(cs), benefitmetrics(bhat, S1, S0))
}

sim_one <- function(dataset, hazard, beta_trt, n, cens, stress, cens_b) {
  ref <- REFS[[dataset]]; times <- TIMES[[cell_key(dataset, hazard, beta_trt, stress, cens)]]
  hp <- TUNED[[dataset]][[hazard]][[as.character(n)]]
  train <- simulate_stress(ref, n, beta_trt, hazard, stress, cens_b); attr(train, "truth") <- NULL
  test <- simulate_stress(ref, N_TEST, beta_trt, hazard, stress, cens_b)
  tr <- attr(test, "truth"); attr(test, "truth") <- NULL
  S_true <- stress_true_surv(tr$X, tr$trt, ref, beta_trt, hazard, stress, times)
  S1 <- stress_true_surv(tr$X, rep(1L, N_TEST), ref, beta_trt, hazard, stress, times[5])[, 1]
  S0 <- stress_true_surv(tr$X, rep(0L, N_TEST), ref, beta_trt, hazard, stress, times[5])[, 1]
  y <- Surv(test$time, test$status)
  d1 <- test; d1$trt <- 1L; d0 <- test; d0$trt <- 0L
  do.call(rbind, lapply(LEARNERS, function(l) {
    t0 <- proc.time()[["elapsed"]]
    v <- tryCatch({
      f <- fit_any(l, train, ref, hp[[l]])
      S <- predict_any(f, test, times)
      bhat <- predict_any(f, d1, times[5])[, 1] - predict_any(f, d0, times[5])[, 1]
      if (anyNA(S) || anyNA(bhat)) stop("NA predictions")
      score_all(S, S_true, y, times, bhat, S1, S0)
    }, error = function(e) NA_ROW)
    v["secs"] <- proc.time()[["elapsed"]] - t0
    data.frame(learner = l, metric = names(v), value = unname(v))
  }))
}

if (sys.nframe() == 0L) {
  a <- commandArgs(TRUE)
  reps <- if (length(a) >= 1) as.integer(a[1]) else 90L
  cores <- if (length(a) >= 2) as.integer(a[2]) else 2L
  ttag <- if (length(a) >= 3) a[3] else "v2"
  TUNED <- readRDS(sprintf("results/tuned_%s.rds", ttag))
  grid <- prepare_grid(GRID); rownames(grid) <- NULL
  TIMES <- grid_times(grid)
  tag <- sprintf("v2_%dreps%s", reps, if (ttag == "v2") "" else paste0("_", ttag))
  saveRDS(list(grid = grid, times = TIMES, learners = LEARNERS, tuned = TUNED, base = BASE,
               versions = pkg_versions()), sprintf("results/design_%s.rds", tag))
  raw <- runsim(n = reps, grid_params = grid, sim_func = sim_one, seed = 2028, num_cores = cores,
                show_progress = FALSE, backend = "fork", on_error = "warn",
                save_path = sprintf("results/checkpoint_%s.rds", tag), checkpoint_every = 5L)
  saveRDS(raw, sprintf("results/sim_raw_%s.rds", tag))
  write.csv(raw, sprintf("results/sim_raw_%s.csv", tag), row.names = FALSE)
}

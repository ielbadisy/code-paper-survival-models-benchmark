# armA-rct/dca-run.R -- treatment decision curves against the truth: rule "treat if predicted benefit > delta",
# delta = harm or cost of treatment in survival points at t*. Net value V(d) - delta * P(d = 1); gain over the
# better of treat-all and treat-none; headroom of the optimal rule 1{b > delta} (learner "optimal").
# Subset of the base design: n = 200, 30% censoring, both trials, 4 hazards, 3 log HRs (24 cells).
# In: common.R, results/tuned_v2.rds. Out: results/dca_raw.rds, results/analysis_v2_90reps/dca.csv
# Run from armA-rct/:  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript dca-run.R <reps> <cores>
source("common.R")
DELTA <- seq(0, 0.1, 0.01)
# gain of rule d at each delta, scored with the true benefit b of every test patient
gain <- function(d, b) sapply(seq_along(DELTA), function(k) mean(d[, k] * (b - DELTA[k])) - max(mean(b) - DELTA[k], 0))
row <- function(l, v) data.frame(learner = l, metric = sprintf("gain_%02d", round(100 * DELTA)), value = v)

sim_one <- function(dataset, hazard, beta_trt, n, cens, stress, cens_b) {
  ref <- REFS[[dataset]]; ts <- TIMES[[cell_key(dataset, hazard, beta_trt, stress, cens)]][5]
  hp <- TUNED[[dataset]][[hazard]][[as.character(n)]]
  train <- simulate_stress(ref, n, beta_trt, hazard, stress, cens_b); attr(train, "truth") <- NULL
  test <- simulate_stress(ref, N_TEST, beta_trt, hazard, stress, cens_b)
  X <- attr(test, "truth")$X; attr(test, "truth") <- NULL
  b <- stress_true_surv(X, rep(1L, N_TEST), ref, beta_trt, hazard, stress, ts)[, 1] -
       stress_true_surv(X, rep(0L, N_TEST), ref, beta_trt, hazard, stress, ts)[, 1]
  d1 <- test; d1$trt <- 1L; d0 <- test; d0$trt <- 0L
  rbind(row("optimal", gain(outer(b, DELTA, ">"), b)), do.call(rbind, lapply(LEARNERS, function(l) {
    v <- tryCatch({
      f <- fit_any(l, train, ref, hp[[l]]); bh <- predict_any(f, d1, ts)[, 1] - predict_any(f, d0, ts)[, 1]
      gain(outer(bh, DELTA, ">"), b)
    }, error = function(e) rep(NA_real_, length(DELTA)))
    row(l, v)
  })))
}

a <- commandArgs(TRUE); reps <- as.integer(a[1]); cores <- as.integer(a[2])
TUNED <- readRDS("results/tuned_v2.rds")
grid <- prepare_grid(base_grid[base_grid$n == 200 & base_grid$cens == 0.3, ]); rownames(grid) <- NULL
TIMES <- grid_times(grid)
raw <- runsim(n = reps, grid_params = grid, sim_func = sim_one, seed = 4031, num_cores = cores,
              show_progress = FALSE, backend = "fork", on_error = "warn")
saveRDS(raw, "results/dca_raw.rds")
w <- aggregate(value ~ dataset + hazard + beta_trt + learner + metric, raw, mean)
write.csv(w, "results/analysis_v2_90reps/dca.csv", row.names = FALSE)

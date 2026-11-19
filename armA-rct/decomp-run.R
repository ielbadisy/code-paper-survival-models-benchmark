# armA-rct/decomp-run.R -- empirical check of the benefit error decomposition
#   E(bhat - b)^2 = 2 * PE * (1 - rho),  PE = (E e1^2 + E e0^2) / 2,  rho = E(e1 e0) / PE,
# with e_a = S_hat(t* | x, a) - S(t* | x, a) on the test trial (both arms for every patient).
# Subset of the base design: n = 200, 30% censoring, both trials, 4 hazards, 3 log HRs (24 cells).
# In: common.R, results/tuned_v2.rds. Out: results/decomp_raw.rds, results/analysis_v2_90reps/decomp.csv
# Run from armA-rct/:  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript decomp-run.R <reps> <cores>

source("common.R")

sim_one <- function(dataset, hazard, beta_trt, n, cens, stress, cens_b) {
  ref <- REFS[[dataset]]; ts <- TIMES[[cell_key(dataset, hazard, beta_trt, stress, cens)]][5]
  hp <- TUNED[[dataset]][[hazard]][[as.character(n)]]
  train <- simulate_stress(ref, n, beta_trt, hazard, stress, cens_b); attr(train, "truth") <- NULL
  test <- simulate_stress(ref, N_TEST, beta_trt, hazard, stress, cens_b)
  X <- attr(test, "truth")$X; attr(test, "truth") <- NULL
  S1 <- stress_true_surv(X, rep(1L, N_TEST), ref, beta_trt, hazard, stress, ts)[, 1]
  S0 <- stress_true_surv(X, rep(0L, N_TEST), ref, beta_trt, hazard, stress, ts)[, 1]
  d1 <- test; d1$trt <- 1L; d0 <- test; d0$trt <- 0L
  do.call(rbind, lapply(LEARNERS, function(l) {
    v <- tryCatch({
      f <- fit_any(l, train, ref, hp[[l]])
      e1 <- predict_any(f, d1, ts)[, 1] - S1; e0 <- predict_any(f, d0, ts)[, 1] - S0
      c(m11 = mean(e1^2), m00 = mean(e0^2), m10 = mean(e1 * e0), bmse = mean((e1 - e0)^2))
    }, error = function(e) c(m11 = NA, m00 = NA, m10 = NA, bmse = NA))
    data.frame(learner = l, metric = names(v), value = unname(v))
  }))
}

a <- commandArgs(TRUE); reps <- as.integer(a[1]); cores <- as.integer(a[2])
TUNED <- readRDS("results/tuned_v2.rds")
grid <- prepare_grid(base_grid[base_grid$n == 200 & base_grid$cens == 0.3, ]); rownames(grid) <- NULL
TIMES <- grid_times(grid)
raw <- runsim(n = reps, grid_params = grid, sim_func = sim_one, seed = 3029, num_cores = cores,
              show_progress = FALSE, backend = "fork", on_error = "warn")
saveRDS(raw, "results/decomp_raw.rds")
# moments averaged over test patients and replicates, then rho per cell
w <- aggregate(value ~ dataset + hazard + beta_trt + learner + metric, raw, mean)
w <- reshape(w, idvar = c("dataset", "hazard", "beta_trt", "learner"), timevar = "metric", direction = "wide")
names(w) <- sub("^value\\.", "", names(w))
w <- transform(w, pe = (m11 + m00) / 2, rho = m10 / ((m11 + m00) / 2))
write.csv(w, "results/analysis_v2_90reps/decomp.csv", row.names = FALSE)

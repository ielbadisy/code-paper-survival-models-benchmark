# armB-real/actg175-dca.R -- decision curves on ACTG 175, same 50 splits and frozen tuning as actg175-run.R.
#   treatment: rule "treat if predicted benefit > delta", delta = cost of treatment in survival points at t*;
#              net value p * S_KM(treated, d = 1) + (1 - p) * S_KM(control, d = 0) - delta * p, minus the
#              better of treat-all (S_KM(treated) - delta) and treat-none (S_KM(control)).
#   prognostic: IPCW net benefit of flagging risk at t* above a threshold (survalis decision curve).
# In: results/tuned_actg175.rds. Out: results/raw_actg175_dca.rds
# Run from armB-real/:  OMP_NUM_THREADS=1 Rscript actg175-dca.R <cores>
source("actg175-run.R")
DELTA <- seq(0, 0.2, 0.01); PT <- seq(0.05, 0.5, 0.05)

dca_score <- function(S, bhat, test) {
  y <- Surv(test$time, test$status); s1 <- km(test[test$trt == 1, ]); s0 <- km(test[test$trt == 0, ])
  g <- sapply(DELTA, function(dl) val(test, bhat > dl) - dl * mean(bhat > dl) - max(s1 - dl, s0))
  nb <- survalis:::dca_survmat(y, S, TSTAR, PT)$curve  # not exported; same IPCW net benefit as interpret(method = "dca")
  c(setNames(g, sprintf("gain_%02d", round(100 * DELTA))), setNames(nb$model, sprintf("nb_%02d", round(100 * PT))),
    setNames(nb$treat_all, sprintf("nball_%02d", round(100 * PT))), treat_share = mean(bhat > 0))
}

dca_one <- function(dataset) {
  s <- split3(); d1 <- s$test; d1$trt <- 1L; d0 <- s$test; d0$trt <- 0L
  do.call(rbind, lapply(LEARNERS, function(l) {
    v <- tryCatch({
      f <- fitl(l, s$train, TUNED[[l]]); S <- predict_any(f, s$test, TIMES); colnames(S) <- paste0("t=", TIMES)
      bhat <- predict_any(f, d1, TSTAR)[, 1] - predict_any(f, d0, TSTAR)[, 1]
      dca_score(S, bhat, s$test)
    }, error = function(e) c(gain_00 = NA_real_))
    data.frame(learner = l, metric = names(v), value = unname(v))
  }))
}

cores <- as.integer(commandArgs(TRUE)[1]); TUNED <- readRDS("results/tuned_actg175.rds")
raw <- runsim(n = 50L, grid_params = data.frame(dataset = "actg175"), sim_func = dca_one, seed = 2028,
              num_cores = cores, show_progress = FALSE, backend = "fork", on_error = "warn")
saveRDS(list(raw = raw, delta = DELTA, pt = PT, tstar = TSTAR), "results/raw_actg175_dca.rds")

# armA-rct/example-res.R -- real-data example on the two reference trials.
# (1) repeated 10-fold CV of every Arm A learner (prediction metrics),
# (2) RMST difference via survalis estimate() as an ATE sanity check,
# (3) data-only benefit calibration (KM arm difference by quintile of
#     predicted benefit, out of fold).
# Run from armA-rct/:  Rscript example-res.R <cores>

suppressMessages({ library(survalis); library(survival) })
source("sim-run.R")          # LEARNERS, fit_any, predict_any, REFS
args <- commandArgs(TRUE)
cores <- if (length(args)) as.integer(args[1]) else 2L

real_trial <- function(name) {
  if (name == "pbc") {
    d <- subset(survival::pbc, !is.na(trt))
    d$status <- as.integer(d$status == 2)
    d$trt <- as.integer(d$trt == 1)
  } else {
    d <- read.csv("../data-ref/adv_prostate_ca/adv_prostate_ca.csv")
    d$trt <- as.integer(d$rx == 1)
    d$time <- d$survtime; d$status <- d$cens
  }
  X <- REFS[[name]]$X
  data.frame(time = d$time, status = d$status, trt = d$trt, X)
}

cv_one <- function(d, ref, learner, folds, times) {
  S <- matrix(NA_real_, nrow(d), length(times)); bhat <- rep(NA_real_, nrow(d))
  for (k in sort(unique(folds))) {
    tr <- d[folds != k, ]; te <- d[folds == k, ]
    f <- tryCatch(fit_any(learner, tr, ref), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    S[folds == k, ] <- predict_any(f, te, times)
    d1 <- te; d1$trt <- 1L; d0 <- te; d0$trt <- 0L
    bhat[folds == k] <- predict_any(f, d1, times[2])[, 1] - predict_any(f, d0, times[2])[, 1]
  }
  list(S = S, bhat = bhat)
}

out <- list(); cal <- list(); ate <- list()
for (name in names(REFS)) {
  d <- real_trial(name); ref <- REFS[[name]]
  times <- unname(quantile(d$time[d$status == 1], c(0.25, 0.5, 0.75)))
  y <- Surv(d$time, d$status)
  for (r in 1:5) {
    set.seed(100 + r); folds <- sample(rep_len(1:10, nrow(d)))
    res <- parallel::mclapply(LEARNERS, function(l) cv_one(d, ref, l, folds, times), mc.cores = cores)
    for (j in seq_along(LEARNERS)) if (!is.null(res[[j]])) {
      v <- score_matrix(res[[j]]$S, y, times, 2)
      out[[length(out) + 1]] <- data.frame(dataset = name, rep = r, learner = LEARNERS[j],
                                           metric = names(v), value = unname(v))
      if (r == 1) cal[[length(cal) + 1]] <- cbind(dataset = name, learner = LEARNERS[j],
        observedbenefitcal(res[[j]]$bhat, d$time, d$status, d$trt, times[2]))
    }
  }
  ate[[name]] <- tryCatch(estimate(Surv(time, status) ~ ., d, model = "coxph", estimand = "RMST_diff",
                                   treatment = "trt", tau = times[3], R = 200),
                          error = function(e) conditionMessage(e))
}
dir.create("results", showWarnings = FALSE)
cvres <- do.call(rbind, out)
write.csv(aggregate(value ~ dataset + learner + metric, cvres, mean), "results/example_cv.csv", row.names = FALSE)
write.csv(do.call(rbind, cal), "results/example_benefit_cal.csv", row.names = FALSE)
saveRDS(ate, "results/example_ate.rds")

# armB-real/actg175-run.R -- real-trial illustration: the 14 Arm A learners on ACTG 175
# (Hammer et al. 1996), zidovudine alone (arm 0) vs zidovudine + didanosine (arm 1), n = 1054,
# 15 baseline covariates, outcome = AIDS-defining event / CD4 drop / death (cens), time in years.
# No truth in real data: everything is scored on held-out thirds of the trial.
#   tune: 10 random 2:1 splits, each HP config scored by validation IBS; best mean frozen.
#   eval: 50 new 2:1 splits; prediction metrics, survival calibration at t* (deciles, KM),
#         benefit calibration (quintiles of predicted benefit, KM arm difference),
#         policy value of "treat if predicted benefit > 0" (KM within rule-concordant arms).
# Out: results/tuned_actg175.rds, results/raw_actg175.rds
# Run from armB-real/:  OMP_NUM_THREADS=1 Rscript actg175-run.R tune|eval <cores> [splits]

setwd("../armA-rct"); source("common.R"); HP <- local({ source("tune-run.R", local = TRUE); HP })
setwd("../armB-real")

COVS <- c("age", "wtkg", "hemo", "homo", "drugs", "karnof", "oprior", "z30", "preanti", "race", "gender",
          "str2", "symptom", "cd40", "cd80")
TSTAR <- 2
LEARNERS[LEARNERS == "cox_spec"] <- "cox_int"   # no oracle in real data: Cox with every trt interaction
data(ACTG175, package = "BART")
D <- with(ACTG175[ACTG175$arms %in% 0:1, ], data.frame(time = days / 365.25, status = cens, trt = as.integer(arms == 1),
                                                       ACTG175[ACTG175$arms %in% 0:1, COVS]))
REF <- list(X = D[COVS], inter = setNames(vector("list", length(COVS)), COVS))
tau <- unname(quantile(D$time, 0.9))
TIMES <- sort(unique(c(unname(quantile(D$time[D$status == 1 & D$time <= tau], seq(0.1, 0.9, 0.1))), TSTAR)))
K <- match(TSTAR, TIMES)

split3 <- function() { i <- unlist(lapply(split(seq_len(nrow(D)), D$trt), function(j) sample(j, round(length(j) / 3))))
  list(train = D[-i, ], test = D[i, ]) }
fitl <- function(l, train, hp) fit_any(if (l == "cox_int") "cox_spec" else l, train, REF, hp)
km <- function(d, t = TSTAR) if (!nrow(d)) NA_real_ else
  summary(survfit(Surv(time, status) ~ 1, d), times = t, extend = TRUE)$surv
# value of rule d by KM within rule-concordant arms; an empty group (weight p or 1 - p near 0) uses its whole arm
val <- function(test, d) { p <- mean(d); a1 <- d & test$trt == 1; a0 <- !d & test$trt == 0
  p * km(test[if (any(a1)) a1 else test$trt == 1, ]) + (1 - p) * km(test[if (any(a0)) a0 else test$trt == 0, ]) }
bins <- function(v, g) cut(v, unique(quantile(v, seq(0, 1, length.out = g + 1))), include.lowest = TRUE, labels = FALSE)

tune_one <- function(dataset) {
  s <- split3(); y <- Surv(s$test$time, s$test$status)
  do.call(rbind, lapply(names(HP), function(l) do.call(rbind, lapply(seq_along(HP[[l]]), function(k)
    data.frame(learner = l, config = k, ibs = tryCatch({
      S <- predict_any(fitl(l, s$train, HP[[l]][[k]]), s$test, TIMES)
      if (anyNA(S)) NA_real_ else evaluate(S, y = y, times = TIMES, metrics = "ibs")$summary$mean
    }, error = function(e) NA_real_))))))
}

# survival calibration (10 groups of predicted S(t*)), benefit calibration (5 groups of predicted
# benefit), policy value; metric names carry the bin index so runsim can stack them
score_real <- function(S, bhat, test) {
  y <- Surv(test$time, test$status); risk <- 1 - S[, K]
  g <- bins(S[, K], 10); gb <- bins(bhat, 5)
  cs <- sapply(sort(unique(g)), function(k) c(mean(S[g == k, K]), km(test[g == k, ])))
  cb <- sapply(sort(unique(gb)), function(k) { d <- test[gb == k, ]
    c(mean(bhat[gb == k]), km(d[d$trt == 1, ]) - km(d[d$trt == 0, ]), nrow(d)) })
  c(score_matrix(S, y, TIMES, K),
    uno_c = concordance(y ~ risk, reverse = TRUE, timewt = "n/G2", ymax = TSTAR)$concordance,
    dwrange = evaluate(S[, K, drop = FALSE], y = y, times = TSTAR, metrics = DW_RANGE)$summary$mean,
    cal_slope = unname(coef(lm(cs[2, ] ~ cs[1, ]))[2]),
    benefit_cal_slope = unname(coef(lm(cb[2, ] ~ cb[1, ], weights = cb[3, ]))[2]),
    benefit_spread = cb[2, ncol(cb)] - cb[2, 1], treat_share = mean(bhat > 0),
    policy_gain = val(test, bhat > 0) - max(km(test[test$trt == 1, ]), km(test[test$trt == 0, ])),
    setNames(cs[1, ], paste0("calS_pred_", seq_len(ncol(cs)))), setNames(cs[2, ], paste0("calS_obs_", seq_len(ncol(cs)))),
    setNames(cb[1, ], paste0("calB_pred_", seq_len(ncol(cb)))), setNames(cb[2, ], paste0("calB_obs_", seq_len(ncol(cb)))))
}

eval_one <- function(dataset) {
  s <- split3(); d1 <- s$test; d1$trt <- 1L; d0 <- s$test; d0$trt <- 0L
  do.call(rbind, lapply(LEARNERS, function(l) {
    t0 <- proc.time()[["elapsed"]]
    v <- tryCatch({
      f <- fitl(l, s$train, TUNED[[l]])
      S <- predict_any(f, s$test, TIMES)
      bhat <- predict_any(f, d1, TSTAR)[, 1] - predict_any(f, d0, TSTAR)[, 1]
      if (anyNA(S) || anyNA(bhat)) stop("NA predictions")
      score_real(S, bhat, s$test)
    }, error = function(e) c(ibs = NA_real_))
    v["secs"] <- proc.time()[["elapsed"]] - t0
    data.frame(learner = l, metric = names(v), value = unname(v))
  }))
}

if (sys.nframe() == 0L) {
  a <- commandArgs(TRUE); mode <- a[1]; cores <- as.integer(a[2])
  dir.create("results", showWarnings = FALSE); g <- data.frame(dataset = "actg175")
  if (mode == "tune") {
    raw <- runsim(n = 10L, grid_params = g, sim_func = tune_one, seed = 1027, num_cores = cores,
                  show_progress = FALSE, backend = "fork", on_error = "warn")
    agg <- aggregate(ibs ~ learner + config, raw, mean)
    best <- do.call(rbind, lapply(split(agg, agg$learner), function(d) d[which.min(d$ibs), ]))
    TUNED <- setNames(Map(function(l, k) HP[[l]][[k]], best$learner, best$config), best$learner)
    saveRDS(list(raw = raw, best = best, tuned = TUNED, versions = pkg_versions()), "results/tune_raw_actg175.rds")
    saveRDS(TUNED, "results/tuned_actg175.rds")
  } else {
    TUNED <- readRDS("results/tuned_actg175.rds"); n <- if (length(a) >= 3) as.integer(a[3]) else 50L
    raw <- runsim(n = n, grid_params = g, sim_func = eval_one, seed = 2028, num_cores = cores,
                  show_progress = FALSE, backend = "fork", on_error = "warn")
    saveRDS(list(raw = raw, times = TIMES, tstar = TSTAR, tuned = TUNED, versions = pkg_versions()),
            sprintf("results/raw_actg175_%dsplits.rds", n))
  }
}

# R/sim-core.R -- fit/predict/score shared by every sim-run.R in the project.
# Ported to the survalis 2.4.2 API (2026-11-15):
#   fit(formula, data, model = <learner>, spec = list(...), seed = )
#   predict(fit, newdata, times) -> survmat (n x length(times))
#   evaluate(<matrix>, y = Surv(...), times, metrics) -> $summary$mean
# Metric names follow listmetrics(). cindex/auc are taken at max(times);
# brier/ece/survacc need a single time and are scored at t_single.

suppressMessages({ library(survalis); library(survival) })

metric_names <- c("cindex", "auc", "ibs", "iae", "ise", "brier", "ece", "survacc")
truth_metric_names <- c("uno_c", "mae_true", "mae_true_tstar")

fit_learner <- function(learner, formula, train, hyperparams = NULL, seed = NULL) {
  fit(formula, train, model = learner, spec = if (is.null(hyperparams)) list() else hyperparams,
      seed = seed)
}

# metric vector from a prediction matrix S (n x length(times)), single-time
# metrics at column k_single
score_matrix <- function(S, y, times, k_single = ceiling(length(times) / 2)) {
  multi <- evaluate(S, y = y, times = times, metrics = c("cindex", "auc", "ibs", "iae", "ise"))
  single <- evaluate(S[, k_single, drop = FALSE], y = y, times = times[k_single],
                     metrics = c("brier", "ece", "survacc"))
  setNames(c(multi$summary$mean, single$summary$mean), metric_names)
}

# truth-based metrics: Uno C on predicted risk at t_single, and mean absolute
# error of S_hat against the true S(t | x)
score_truth <- function(S, S_true, y, times, k_single = ceiling(length(times) / 2)) {
  risk <- 1 - S[, k_single]
  c(uno_c = concordance(y ~ risk, reverse = TRUE, timewt = "n/G2", ymax = max(times))$concordance,
    mae_true = mean(abs(S - S_true)),
    mae_true_tstar = mean(abs(S[, k_single] - S_true[, k_single])))
}

# backwards-compatible scorer used by paper1-3 (training-quartile time grid)
score_fit <- function(learner, fit, train, test) {
  tryCatch({
    times <- unname(quantile(train$time[train$status == 1], seq(0.2, 0.8, 0.2)))
    S <- as.matrix(predict(fit, newdata = test, times = times))
    v <- score_matrix(S, Surv(test$time, test$status), times, k_single = 3)
    data.frame(metric = names(v), value = unname(v))
  }, error = function(e) data.frame(metric = metric_names, value = NA_real_))
}

fit_predict_score <- function(learner, formula, train, test, hyperparams) {
  out <- tryCatch({
    fit <- fit_learner(learner, formula, train, hyperparams)
    score_fit(learner, fit, train, test)
  }, error = function(e) data.frame(metric = metric_names, value = NA_real_))
  out$learner <- learner
  out
}

# factory for the general factorial arm (paper1/2 and arm B): closes over
# the censoring cache, tuned hyperparameters, and learner set
make_sim_one <- function(censoring_params, tuned_params, learners, n_test = 1000L) {
  function(dgp, censor_rate, n, p, data_type, mechanism) {
    rate0 <- censoring_params$rate0[
      censoring_params$dgp == dgp & censoring_params$censor_rate == censor_rate &
      censoring_params$p == p & censoring_params$data_type == data_type &
      censoring_params$mechanism == mechanism
    ]
    seed <- sample.int(1e6, 1)
    train <- simulate_dataset(dgp, n, p, data_type, censor_rate, mechanism, rate0 = rate0, seed = seed)
    test <- simulate_dataset(dgp, n_test, p, data_type, censor_rate, mechanism, rate0 = rate0, seed = seed + 500000L)
    covs <- setdiff(names(train), c("time", "status"))
    formula <- as.formula(paste("Surv(time, status) ~", paste(covs, collapse = " + ")))
    times <- unname(quantile(train$time[train$status == 1], seq(0.2, 0.8, 0.2)))
    S_true <- true_surv(dgp, test[covs], times)
    y <- Surv(test$time, test$status)
    dgp_params <- tuned_params[[dgp]]
    do.call(rbind, lapply(learners, function(l) {
      v <- tryCatch({
        f <- fit_learner(l, formula, train, dgp_params[[l]])
        S <- as.matrix(predict(f, newdata = test, times = times))
        c(score_matrix(S, y, times, 3), score_truth(S, S_true, y, times, 3))
      }, error = function(e) setNames(rep(NA_real_, 11), c(metric_names, truth_metric_names)))
      data.frame(learner = l, metric = names(v), value = unname(v))
    }))
  }
}

# factory for paper3 (missing data): unchanged design, new API underneath
make_sim_one_missingness <- function(censoring_params, tuned_params, learners, n_test = 1000L,
                                      dgp = "ph_weibull", n = 700L, p = 10L, data_type = "mixed",
                                      censor_rate = 0.4, mechanism = "independent") {
  rate0 <- censoring_params$rate0[
    censoring_params$dgp == dgp & censoring_params$censor_rate == censor_rate &
    censoring_params$p == p & censoring_params$data_type == data_type &
    censoring_params$mechanism == mechanism
  ]
  dgp_params <- tuned_params[[dgp]]

  function(missing_rate, miss_mechanism) {
    seed <- sample.int(1e6, 1)
    train <- simulate_dataset(dgp, n, p, data_type, censor_rate, mechanism, rate0 = rate0, seed = seed)
    test <- simulate_dataset(dgp, n_test, p, data_type, censor_rate, mechanism, rate0 = rate0, seed = seed + 500000L)
    covs <- setdiff(names(train), c("time", "status"))
    formula <- as.formula(paste("Surv(time, status) ~", paste(covs, collapse = " + ")))

    test_amp <- ampute_covariates(test, prop = missing_rate, mechanism = miss_mechanism, seed = seed + 1L)
    cc <- test_amp[stats::complete.cases(test_amp), , drop = FALSE]
    enough_cc <- nrow(cc) >= 20

    out <- do.call(rbind, lapply(learners, function(l) {
      fit <- tryCatch(fit_learner(l, formula, train, dgp_params[[l]]), error = function(e) NULL)
      na <- data.frame(metric = metric_names, value = NA_real_)
      ref_l <- if (is.null(fit)) na else score_fit(l, fit, train, test)
      miss_l <- if (is.null(fit) || !enough_cc) na else score_fit(l, fit, train, cc)
      ref_l$learner <- l; ref_l$arm <- "reference"; ref_l$n_eval <- nrow(test)
      miss_l$learner <- l; miss_l$arm <- "missing"; miss_l$n_eval <- nrow(cc)
      rbind(ref_l, miss_l)
    }))
    out$missing_rate <- missing_rate
    out$miss_mechanism <- miss_mechanism
    out
  }
}

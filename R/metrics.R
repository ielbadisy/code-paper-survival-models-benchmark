# R/metrics.R -- performance and benefit measures (survalis evaluate()).
# Prediction: cindex/auc/ibs/iae/ise over the time grid, brier/ece/survacc at t* (column k).
# Benefit at t*: b(x) = S(t* | x, 1) - S(t* | x, 0), predicted by scoring each patient with trt = 1 and 0.
suppressMessages({ library(survalis); library(survival) })
metric_names <- c("cindex", "auc", "ibs", "iae", "ise", "brier", "ece", "survacc")
truth_metric_names <- c("uno_c", "mae_true", "mae_true_tstar")

score_matrix <- function(S, y, times, k) {
  multi <- evaluate(S, y = y, times = times, metrics = c("cindex", "auc", "ibs", "iae", "ise"))
  single <- evaluate(S[, k, drop = FALSE], y = y, times = times[k], metrics = c("brier", "ece", "survacc"))
  setNames(c(multi$summary$mean, single$summary$mean), metric_names)
}
# Uno C on risk at t*, mean |S_hat - S_true| over the grid and at t*
score_truth <- function(S, S_true, y, times, k) {
  risk <- 1 - S[, k]
  c(uno_c = concordance(y ~ risk, reverse = TRUE, timewt = "n/G2", ymax = max(times))$concordance,
    mae_true = mean(abs(S - S_true)), mae_true_tstar = mean(abs(S[, k] - S_true[, k])))
}
# truth-based c-for-benefit: share of random pairs ordered alike by predicted and true benefit
# (truth ties skipped, prediction ties count 1/2)
.benefitconcordance <- function(bhat, btrue, max_pairs = 2e5) {
  n <- length(bhat); i <- sample.int(n, max_pairs, TRUE); j <- sample.int(n, max_pairs, TRUE)
  keep <- btrue[i] != btrue[j]
  if (!any(keep)) return(NA_real_)
  i <- i[keep]; j <- j[keep]; s <- sign(bhat[i] - bhat[j])
  mean(ifelse(s == 0, 0.5, s == sign(btrue[i] - btrue[j])))
}
# bhat: predicted benefit; S1, S0: true S(t* | x, 1 / 0). Rule: treat if bhat > delta, valued on the truth.
benefitmetrics <- function(bhat, S1, S0, delta = 0) {
  btrue <- S1 - S0; val <- function(treat) mean(ifelse(treat, S1, S0))
  value <- val(bhat > delta)
  fitc <- if (var(bhat) > 0) coef(lm(btrue ~ bhat)) else c(NA, NA)
  c(benefit_rmse = sqrt(mean((bhat - btrue)^2)), benefit_bias = mean(bhat - btrue),
    benefit_cal_int = unname(fitc[1]), benefit_cal_slope = unname(fitc[2]),
    benefit_cindex = .benefitconcordance(bhat, btrue),
    policy_regret = val(btrue > delta) - value,           # survival lost vs the oracle rule
    policy_gain = value - max(mean(S1), mean(S0)))        # gain over treat-all / treat-none
}

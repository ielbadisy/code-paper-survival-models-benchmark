# R/benefit.R -- individual treatment-benefit metrics for simulated trials.
#
# Benefit at horizon t*: b(x) = S(t* | x, trt = 1) - S(t* | x, trt = 0),
# the absolute gain in survival probability from treatment. A learner's
# predicted benefit is obtained by predicting each test patient twice, with
# trt set to 1 and to 0 (counterfactual prediction). The true benefit comes
# from rct_true_surv() in R/dgp-rct.R. Positive benefit = treatment helps.

# predicted benefit at times `times` for a fitted survalis model
predictbenefit <- function(fit, newdata, times) {
  d1 <- newdata; d1$trt <- 1L
  d0 <- newdata; d0$trt <- 0L
  as.matrix(predict(fit, newdata = d1, times = times)) -
    as.matrix(predict(fit, newdata = d0, times = times))
}

# true benefit (n x length(times)) from the DGP
truebenefit <- function(X, ref, beta_trt, hazard, times, centre = TRUE) {
  n <- nrow(X)
  rct_true_surv(X, rep(1L, n), ref, beta_trt, hazard, times, centre) -
    rct_true_surv(X, rep(0L, n), ref, beta_trt, hazard, times, centre)
}

# probability that two random patients are ordered the same way by
# predicted and true benefit (ties in the truth skipped, ties in the
# prediction count 1/2); a truth-based c-for-benefit
.benefitconcordance <- function(bhat, btrue, max_pairs = 2e5) {
  n <- length(bhat)
  i <- sample.int(n, max_pairs, TRUE); j <- sample.int(n, max_pairs, TRUE)
  keep <- btrue[i] != btrue[j]
  if (!any(keep)) return(NA_real_)
  i <- i[keep]; j <- j[keep]
  s_true <- sign(btrue[i] - btrue[j]); s_hat <- sign(bhat[i] - bhat[j])
  mean(ifelse(s_hat == 0, 0.5, s_hat == s_true))
}

# value of the rule "treat if predicted benefit > delta", on the truth scale:
# mean true survival at t* when each patient gets the treatment the rule picks
.policyvalue <- function(treat, S1, S0) mean(ifelse(treat, S1, S0))

# all benefit metrics at a single horizon t* (vectors of length n)
# bhat: predicted benefit; S1, S0: true S(t*|x, trt = 1 / 0)
benefitmetrics <- function(bhat, S1, S0, delta = 0) {
  btrue <- S1 - S0
  oracle <- .policyvalue(btrue > delta, S1, S0)
  best_single <- max(mean(S1), mean(S0))        # better of treat-all / treat-none
  value <- .policyvalue(bhat > delta, S1, S0)
  fitc <- if (stats::var(bhat) > 0) stats::coef(stats::lm(btrue ~ bhat)) else c(NA, NA)
  c(benefit_rmse = sqrt(mean((bhat - btrue)^2)),
    benefit_bias = mean(bhat - btrue),
    benefit_cal_int = unname(fitc[1]),           # truth-based calibration intercept
    benefit_cal_slope = unname(fitc[2]),         # ... and slope (1 = perfect)
    benefit_cindex = .benefitconcordance(bhat, btrue),
    policy_regret = oracle - value,              # survival lost vs the oracle rule
    policy_gain = value - best_single)           # gain over treating everyone or no one
}

# data-only benefit calibration (what an analyst could compute without the
# truth): Kaplan-Meier arm difference at t* within groups of predicted benefit
observedbenefitcal <- function(bhat, time, status, trt, tstar, groups = 5) {
  g <- cut(bhat, unique(stats::quantile(bhat, seq(0, 1, length.out = groups + 1))),
           include.lowest = TRUE, labels = FALSE)
  do.call(rbind, lapply(sort(unique(g)), function(k) {
    i <- g == k
    km <- function(a) {
      ii <- i & trt == a
      if (sum(ii) < 5) return(NA_real_)
      s <- summary(survival::survfit(survival::Surv(time[ii], status[ii]) ~ 1),
                   times = tstar, extend = TRUE)$surv
      if (length(s)) s else NA_real_
    }
    data.frame(group = k, predicted = mean(bhat[i]), observed = km(1) - km(0), n = sum(i))
  }))
}

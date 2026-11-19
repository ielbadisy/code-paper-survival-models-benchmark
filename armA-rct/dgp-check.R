# dgp-check.R -- large-sample checks of R/dgp.R before the main run (paper Section 2.5).
# 16 cells: 2 trials x 4 stressors x 2 hazards, beta_trt = -0.4, 30% censoring, n = 20 000.
# Output: results/dgp_check.csv (censored share, max |truth - latent empirical S|, max |KM - truth|
# at the event-time quartiles). Run from armA-rct/: Rscript dgp-check.R
source("../R/dgp.R")
set.seed(11)
chk <- function(ds, st, hz) {
  ref <- load_reference(ds)
  b <- calibrate_stress_censoring(ref, -0.4, hz, st, 0.3, n_pilot = 5e4)
  d <- simulate_stress(ref, 20000, -0.4, hz, st, b); tr <- attr(d, "truth")
  tt <- quantile(d$time[d$status == 1], c(.25, .5, .75))
  true <- colMeans(stress_true_surv(tr$X, tr$trt, ref, -0.4, hz, st, tt))
  # KM is biased under covariate-dependent censoring, so the truth is also checked against latent times
  lat <- stress_draw_times(tr$X, tr$trt, ref, -0.4, hz, st)
  km <- summary(survfit(Surv(time, status) ~ 1, d), times = tt)$surv
  data.frame(dataset = ds, stress = st, hazard = hz, cens = 1 - mean(d$status),
             truth_vs_latent = max(abs(sapply(tt, function(t) mean(lat > t)) - true)),
             km_vs_truth = max(abs(km - true)))
}
g <- expand.grid(hazard = c("const", "nph"), stress = c("none", "nle", "tve", "censdep"),
                 dataset = c("pbc", "prostate"), stringsAsFactors = FALSE)
write.csv(do.call(rbind, Map(chk, g$dataset, g$stress, g$hazard)), "results/dgp_check.csv", row.names = FALSE)

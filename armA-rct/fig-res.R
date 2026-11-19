# armA-rct/fig-res.R -- paper figures and tables for the Arm A v2 run (90 eval reps).
# Main text uses the core design (n = 200, 30% censoring, 24 base scenarios); Fig 1 adds n, Fig 4 the
# stressors against their matched base scenarios. Full-design tables (*_full.csv) go to S1.
# In: results/sim_raw_<tag>.rds, results/analysis_<tag>/*.csv (from explore-res.R)
# Out: results/figures_<tag>/fig{1..5}.{pdf,png}, tab1_core.csv, tab1_full.csv, tab2_core.csv, tab2_full.csv
# Run from armA-rct/:  Rscript fig-res.R [tag]

source("../R/res-utils.R")
a <- commandArgs(TRUE); TAG <- if (length(a)) a[1] else "v2_90reps"
IN <- sprintf("results/analysis_%s", TAG); OUT <- sprintf("results/figures_%s", TAG)
dir.create(OUT, showWarnings = FALSE)
rd <- function(f) read.csv(file.path(IN, f), stringsAsFactors = FALSE)


r <- readRDS(sprintf("results/sim_raw_%s.rds", TAG))
cell <- c("dataset", "hazard", "beta_trt", "n", "cens", "stress")
W <- widen(r, c(cell, "ID", "learner"))
Z95 <- qnorm(0.975)  # 95% Monte Carlo intervals: estimate +/- Z95 * MCSE
is_core <- function(d) d$stress == "none" & d$n == 200 & d$cens == 0.3
# paired difference from Cox PH within replicate: mean and MCSE by `by`
cx <- W[W$learner == "coxph", ]; ix <- match(do.call(paste, W[c(cell, "ID")]), do.call(paste, cx[c(cell, "ID")]))
pdiff <- function(keep, by, m) {
  D <- cbind(W[keep, c(by, "learner")], d = W[keep, m] - cx[ix[keep], m]); D <- D[D$learner != "coxph", ]
  a <- aggregate(d ~ ., D, mean); s <- aggregate(d ~ ., D, function(v) sd(v) / sqrt(length(v)))
  setNames(cbind(a, s$d), c(by, "learner", m, paste0(m, "_mcse")))
}

# ---- Table 1: overall performance (mean over cells; slopes as median) ------
M1 <- c("mae_true", "ibs", "cindex", "dwrange", "benefit_rmse", "benefit_cindex", "policy_regret", "secs")
# MCSE of a mean over a fixed set of scenarios: scenarios are design points, not a sample
mcse <- function(v, g) sqrt(sum(tapply(v, g, function(x) var(x, na.rm = TRUE) / sum(!is.na(x))))) /
  length(unique(g))
tab1 <- function(V) {
  cm <- aggregate(V[M1], V[c(cell, "learner")], mean, na.rm = TRUE)
  g <- do.call(paste, V[cell]); ii <- split(seq_len(nrow(V)), V$learner)
  se <- function(m) sapply(ii, function(i) mcse(V[i, m], g[i]))
  t1 <- aggregate(cm[M1], cm["learner"], mean, na.rm = TRUE)
  t1$cal_slope_true <- tapply(V$cal_slope_true, V$learner, median, na.rm = TRUE)[t1$learner]
  t1$benefit_cal_slope <- tapply(V$benefit_cal_slope, V$learner, median, na.rm = TRUE)[t1$learner]
  t1$rank_pred <- rank(t1$mae_true); t1$rank_benefit <- rank(t1$benefit_rmse)
  t1$mae_true_mcse <- signif(se("mae_true")[t1$learner], 2)
  t1$benefit_rmse_mcse <- signif(se("benefit_rmse")[t1$learner], 2)
  t1 <- t1[order(t1$mae_true), ]
  t1 <- cbind(model = lab(t1$learner), family = as.character(fam(t1$learner)), t1)
  num <- sapply(t1, is.numeric) & !grepl("mcse", names(t1)); t1[num] <- lapply(t1[num], round, 3); t1
}
write.csv(tab1(W), file.path(OUT, "tab1_full.csv"), row.names = FALSE)
t1 <- tab1(W[is_core(W), ]); write.csv(t1, file.path(OUT, "tab1_core.csv"), row.names = FALSE)

# ---- Fig 1: accuracy vs truth relative to Cox PH, by training size ---------
p1d <- pdiff(W$stress == "none" & W$cens == 0.3, "n", "mae_true")
write.csv(p1d, file.path(OUT, "fig1_data.csv"), row.names = FALSE)
p1d$model <- factor(lab(p1d$learner), lab(t1$learner)); p1d$family <- fam(p1d$learner)
f1 <- ggplot(p1d, aes(factor(n), mae_true, group = model, colour = family)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  geom_errorbar(aes(ymin = mae_true - Z95 * mae_true_mcse, ymax = mae_true + Z95 * mae_true_mcse), width = 0.15) +
  facet_wrap(~ model, nrow = 3) + scale_colour_manual(values = FAM_COL, name = NULL) +
  labs(x = "Training trial size n", y = "Difference in mean |S_hat - S_true| vs Cox PH\n(below 0 = more accurate)",
       title = "Flexible learners help only in small trials; their edge vanishes or reverses by n = 400") + th
save(f1, "fig1_accuracy_by_n", 9, 6.5)

# ---- Fig 2: learner x metric biplot, rank agreement per cell ----------------
# metrics sign-flipped so every arrow points to better performance
SIGN <- c(mae_true = -1, ibs = -1, cindex = 1, uno_c = 1, dwrange = -1, benefit_rmse = -1,
          benefit_cindex = 1, policy_regret = -1)
ML <- c(mae_true = "Low truth error", ibs = "Low IBS", cindex = "High Harrell C", uno_c = "High Uno C",
        dwrange = "Low decision loss", benefit_rmse = "Low benefit RMSE", benefit_cindex = "High c-for-benefit",
        policy_regret = "Low policy regret")
bm <- aggregate(W[is_core(W), names(SIGN)], W[is_core(W), c(cell, "learner")], mean, na.rm = TRUE)
bm <- aggregate(bm[names(SIGN)], bm["learner"], mean, na.rm = TRUE)
Z <- scale(sweep(as.matrix(bm[names(SIGN)]), 2, SIGN, `*`)); rownames(Z) <- bm$learner
pc <- prcomp(Z); ve <- pc$sdev^2 / sum(pc$sdev^2)
write.csv(data.frame(axis = seq_along(ve), var_explained = round(ve, 3)), file.path(OUT, "biplot_variance.csv"),
          row.names = FALSE)
sc <- data.frame(pc$x[, 1:2], learner = rownames(Z))
sc$model <- ifelse(sc$learner == "cox_spec", "Cox, true form (oracle)", lab(sc$learner))
sc$family <- fam(sc$learner)
ld <- data.frame(pc$rotation[, 1:2] * max(abs(sc[1:2])) / max(abs(pc$rotation[, 1:2])) * 0.9,
                 metric = ML[rownames(pc$rotation)])
ld$type <- ifelse(grepl("^benefit|policy", rownames(pc$rotation)), "Benefit", "Prediction")
tx <- rbind(data.frame(x = sc$PC1, y = sc$PC2, label = sc$model, col = FAM_COL[as.integer(sc$family)]),
            data.frame(x = ld$PC1 * 1.05, y = ld$PC2 * 1.05, label = ld$metric, col = "grey20"))
f2a <- ggplot() + geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_segment(data = ld, aes(0, 0, xend = PC1, yend = PC2, linetype = type), colour = "grey35",
               arrow = arrow(length = unit(5, "pt"))) +
  geom_point(data = sc, aes(PC1, PC2, colour = family), size = 2.4) +
  # one repel layer for learner and measure labels so they avoid each other
  ggrepel::geom_text_repel(data = tx, aes(x, y, label = label), colour = tx$col, size = 2.8, seed = 1,
                           box.padding = 0.35, min.segment.length = 0.3, segment.colour = "grey70",
                           max.overlaps = Inf) +
  scale_colour_manual(values = FAM_COL, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.05))) +
  scale_linetype_manual(values = c(Benefit = 1, Prediction = 2), name = NULL) + coord_equal() +
  labs(x = sprintf("Axis 1 (%.0f%%)", 100 * ve[1]), y = sprintf("Axis 2 (%.0f%%)", 100 * ve[2]),
       title = "A. Prediction and benefit measures point in different directions") + th
pb <- rd("pred_vs_benefit_rank_cor.csv"); pb <- pb[is_core(pb) | pb$stress != "none", ]
write.csv(aggregate(rho_rmse ~ beta_trt + stress, pb, median), file.path(OUT, "fig2b_medians.csv"),
          row.names = FALSE)
SC <- c("Log HR 0", "Log HR -0.4", "Log HR 0.8", "Non-linear", "Time-varying", "Dep. censoring")
pb$scenario <- factor(with(pb, ifelse(stress != "none", c(nle = SC[4], tve = SC[5], censdep = SC[6])[stress],
                                      paste("Log HR", beta_trt))), SC)
f2b <- ggplot(pb, aes(scenario, rho_rmse)) + geom_hline(yintercept = 0, linewidth = 0.3) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = "grey40") +
  geom_point(aes(colour = dataset), position = position_jitter(width = 0.12, seed = 1), size = 1.8) +
  scale_colour_manual(values = c(pbc = "#1B9E77", prostate = "#7570B3"), labels = c("PBC", "Prostate"),
                      name = NULL) +
  labs(x = NULL, y = "Spearman rho across the 14 learners\n(prediction rank vs benefit rank, per scenario)",
       title = "B. Agreement is weakest when the average effect is null") +
  th + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save(f2a + f2b + plot_layout(widths = c(1.5, 1)), "fig2_prediction_vs_benefit", 12, 6.5)

# ---- Fig 3: model selection by data metric ---------------------------------
sr <- rd("selection_raw.csv")
sc <- sr[is_core(sr), ]; gc <- do.call(paste, sc[cell])
s3 <- aggregate(cbind(rho, mae_loss, ben_loss) ~ by + hazard, sc, mean)
s3m <- function(v) mapply(function(b, h) { i <- sc$by == b & sc$hazard == h; mcse(sc[[v]][i], gc[i]) },
                          s3$by, s3$hazard)
s3$rho_mcse <- s3m("rho"); s3$mae_loss_mcse <- s3m("mae_loss")
s3$metric <- factor(c(ibs = "IBS", harrell_c = "Harrell C", uno_c = "Uno C",
                      dwrange = "Decision-weighted loss")[s3$by],
                    c("IBS", "Decision-weighted loss", "Uno C", "Harrell C"))
s3$hazard <- factor(c(dec = "Decreasing", const = "Constant", inc = "Increasing",
                      nph = "Non-PH (crossing)")[s3$hazard],
                    c("Decreasing", "Constant", "Increasing", "Non-PH (crossing)"))
MCOL <- c("IBS" = "#1B7837", "Decision-weighted loss" = "#5AAE61", "Uno C" = "#C2A5CF", "Harrell C" = "#762A83")
dg <- position_dodge(width = 0.3)
g3 <- function(y, yl, tt) ggplot(s3, aes(hazard, .data[[y]], colour = metric, group = metric)) +
  geom_line(linewidth = 0.7, position = dg) + geom_point(size = 2.4, position = dg) +
  geom_errorbar(aes(ymin = .data[[y]] - Z95 * .data[[paste0(y, "_mcse")]],
                    ymax = .data[[y]] + Z95 * .data[[paste0(y, "_mcse")]]), width = 0.2, position = dg) +
  scale_colour_manual(values = MCOL, name = "Selection metric") +
  labs(x = "Baseline hazard shape", y = yl, title = tt) + th
f3 <- g3("rho", "Spearman rho with true accuracy", "A. Does the metric rank learners like the truth?") +
  g3("mae_loss", "Excess error of the picked learner\nvs the truly best one",
     "B. Cost of picking by this metric") +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save(f3, "fig3_model_selection", 10, 4.4)

# ---- Fig 4: robustness under stress (benefit RMSE vs Cox PH) ---------------
# base = the 4 scenarios the stressors are layered on (n 200, 30% cens, beta -0.4, const/crossing)
p4 <- pdiff(W$n == 200 & W$cens == 0.3 & W$beta_trt == -0.4 & W$hazard %in% c("const", "nph"), "stress",
            "benefit_rmse")
write.csv(p4, file.path(OUT, "fig4_data.csv"), row.names = FALSE)
BL <- "Matched base scenarios"
p4$stress <- factor(c(none = BL, nle = "Non-linear effects", tve = "Time-varying effect",
                      censdep = "Dependent censoring")[p4$stress],
                    c(BL, "Non-linear effects", "Time-varying effect", "Dependent censoring"))
o4 <- with(p4[p4$stress == BL, ], lab(learner)[order(benefit_rmse)])
p4$model <- factor(lab(p4$learner), rev(o4)); p4$family <- fam(p4$learner)
f4 <- ggplot(p4, aes(benefit_rmse, model, colour = family)) + geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_errorbar(aes(xmin = benefit_rmse - Z95 * benefit_rmse_mcse, xmax = benefit_rmse + Z95 * benefit_rmse_mcse),
                width = 0, orientation = "y") + geom_point(size = 2) + facet_wrap(~ stress, nrow = 1) +
  scale_colour_manual(values = FAM_COL, name = NULL) +
  labs(x = "Difference in benefit RMSE vs Cox PH (left = better)", y = NULL,
       title = "Benefit gains come from modelling the treatment effect, not from generic flexibility") + th
save(f4, "fig4_stress_benefit", 11, 4.8)

# ---- Fig 5: treatment rules -------------------------------------------------
# prostate core scenarios only: in PBC the headroom is zero (Appendix A)
pr <- is_core(W) & W$dataset == "prostate"
p5 <- aggregate(cbind(policy_regret, policy_gain) ~ learner + beta_trt, W[pr, ], mean)
W5 <- W[pr, ]; g5 <- do.call(paste, W5[cell])
p5$gain_mcse <- mapply(function(l, b) { i <- W5$learner == l & W5$beta_trt == b & !is.na(W5$policy_gain)
  mcse(W5$policy_gain[i], g5[i]) }, p5$learner, p5$beta_trt)
hr <- rd("true_benefit_spread.csv")
hr <- aggregate(headroom ~ beta_trt, hr[hr$dataset == "prostate" & hr$cens == 0.3, ], mean)
write.csv(merge(p5, hr), file.path(OUT, "fig5_data.csv"), row.names = FALSE)
p5$model <- factor(lab(p5$learner), lab(t1$learner[order(-t1$policy_regret)])); p5$family <- fam(p5$learner)
EF <- function(b) factor(paste("Log HR =", b), paste("Log HR =", c(0, -0.4, 0.8)))
p5$effect <- EF(p5$beta_trt); hr$effect <- EF(hr$beta_trt)
f5 <- ggplot(p5, aes(policy_gain * 100, model, colour = family)) + geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_vline(data = hr, aes(xintercept = headroom * 100), linetype = 2, colour = "grey40") +
  geom_errorbar(aes(xmin = 100 * (policy_gain - Z95 * gain_mcse), xmax = 100 * (policy_gain + Z95 * gain_mcse)),
                width = 0, orientation = "y") +
  geom_point(size = 2) + facet_wrap(~ effect, nrow = 1) + scale_colour_manual(values = FAM_COL, name = NULL) +
  labs(x = paste("Gain of the rule over the better of treat-all and treat-none (points of survival at t*)",
                 "Dashed line: headroom, the largest gain any rule can have", sep = "\n"),
       y = NULL, title = "Treatment rules in the prostate trial: gains stay far below the headroom") + th
save(f5, "fig5_policy", 10, 4.6)

# ---- Table 2: selection metrics, core (+ crossing hazards) and full design ----
tab2 <- function(S) {
  t <- aggregate(cbind(rho, mae_loss, ben_loss) ~ by, S, mean); g <- do.call(paste, S[cell])
  m <- function(v) sapply(t$by, function(b) mcse(S[[v]][S$by == b], g[S$by == b]))
  transform(t, rho_mcse = m("rho"), loss_mcse = m("mae_loss"))
}
write.csv(tab2(sr), file.path(OUT, "tab2_full.csv"), row.names = FALSE)
t2 <- tab2(sr[is_core(sr), ]); t2$metric <- s3$metric[match(t2$by, s3$by)]
nph <- s3[s3$hazard == "Non-PH (crossing)", ]; j <- match(t2$by, nph$by)
t2$rho_nph <- nph$rho[j]; t2$loss_nph <- nph$mae_loss[j]
write.csv(t2[order(t2$mae_loss), ], file.path(OUT, "tab2_core.csv"), row.names = FALSE)

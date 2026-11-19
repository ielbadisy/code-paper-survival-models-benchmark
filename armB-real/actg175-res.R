# armB-real/actg175-res.R -- tables and figures for the ACTG 175 illustration.
# In: results/raw_actg175_<n>splits.rds (actg175-run.R eval)
# Out: results/figures/tab_actg175.csv, fig7_actg175_metrics, fig8_actg175_calibration,
#      fig9_actg175_benefit_calibration (.pdf/.png)
# Run from armB-real/:  Rscript actg175-res.R [splits]

suppressMessages({ library(ggplot2); library(patchwork) })
a <- commandArgs(TRUE); N <- if (length(a)) a[1] else "50"
R <- readRDS(sprintf("results/raw_actg175_%ssplits.rds", N)); r <- R$raw
OUT <- "results/figures"; dir.create(OUT, showWarnings = FALSE)

COX <- c("coxph", "cox_int", "glmnet", "coxboost", "rpsurv", "rpsurv_tve")
LAB <- c(coxph = "Cox PH", cox_int = "Cox, all trt interactions", glmnet = "Ridge/lasso Cox", coxboost = "CoxBoost",
         rpsurv = "Royston-Parmar", rpsurv_tve = "Royston-Parmar TVE", nbsurv = "Naive Bayes",
         fastgbm = "GBM (fastgbm)", densemlp = "Neural net (MLP)", bart = "BART", xgboost = "XGBoost",
         rsf = "Random survival forest", orsf = "Oblique RSF", blackboost = "Tree boosting (mboost)")
FAM_COL <- c("Regression" = "#2166AC", "Machine learning" = "#D6604D")
fam <- function(l) factor(ifelse(l %in% COX, "Regression", "Machine learning"), names(FAM_COL))
lab <- function(l) unname(LAB[l])
th <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 11), plot.title.position = "plot")
save <- function(p, f, w, h) {
  ggsave(file.path(OUT, paste0(f, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT, paste0(f, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
}

# ---- table: mean (MCSE) over splits; calibration slopes as medians ---------
s <- r[!grepl("^cal[SB]_", r$metric), ]
W <- reshape(s[c("ID", "learner", "metric", "value")], idvar = c("ID", "learner"), timevar = "metric", direction = "wide")
names(W) <- sub("^value\\.", "", names(W))
M <- c("ibs", "brier", "cindex", "uno_c", "dwrange", "benefit_spread", "treat_share", "policy_gain", "secs")
mn <- aggregate(W[M], W["learner"], mean, na.rm = TRUE)
se <- aggregate(W[M], W["learner"], function(v) sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v))))
names(se)[-1] <- paste0(M, "_mcse")
tb <- merge(mn, se)
tb$cal_slope <- tapply(W$cal_slope, W$learner, median, na.rm = TRUE)[tb$learner]
tb$benefit_cal_slope <- tapply(W$benefit_cal_slope, W$learner, median, na.rm = TRUE)[tb$learner]
tb$fail <- tapply(is.na(W$ibs), W$learner, mean)[tb$learner]
tb <- tb[order(tb$ibs), ]; tb <- cbind(model = lab(tb$learner), tb)
num <- sapply(tb, is.numeric); out <- tb; out[num] <- lapply(out[num], round, 4)
write.csv(out, file.path(OUT, "tab_actg175.csv"), row.names = FALSE)

# ---- fig 7: do data metrics agree, and does prediction track benefit? -----
# benefit by observed spread, not |slope - 1|: shrunk predictions give slope ~1 without accuracy
rk <- within(tb, { r_ibs <- rank(ibs); r_c <- rank(-cindex); r_ben <- rank(-benefit_spread) })
rk$model <- lab(rk$learner); rk$family <- fam(rk$learner)
sc <- function(x, y, xl, yl, tt) ggplot(rk, aes(.data[[x]], .data[[y]], colour = family)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") + geom_point(size = 2.6) +
  ggrepel::geom_text_repel(aes(label = model), size = 2.9, show.legend = FALSE, seed = 1) +
  scale_colour_manual(values = FAM_COL, name = NULL) + coord_equal(xlim = c(1, 14), ylim = c(1, 14)) +
  labs(x = xl, y = yl, title = tt) + th
rho <- function(x, y) sprintf("Spearman rho = %.2f", cor(rk[[x]], rk[[y]], method = "spearman"))
f7 <- sc("r_ibs", "r_c", "Rank by IBS (1 = best)", "Rank by Harrell C", paste("A. IBS vs Harrell C.", rho("r_ibs", "r_c"))) +
  sc("r_ibs", "r_ben", "Rank by IBS (1 = best)", "Rank by observed benefit spread\n(top minus bottom quintile)",
     paste("B. Prediction vs benefit.", rho("r_ibs", "r_ben"))) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save(f7 + plot_annotation(title = sprintf("ACTG 175 (n = 1054), %s random 2:1 splits", N),
                          theme = theme(plot.title = element_text(face = "bold"))), "fig7_actg175_metrics", 11, 5.8)

# ---- calibration curves: bin means averaged over splits, 10-90% band -------
cal <- function(px) {
  c <- r[grepl(paste0("^", px, "_"), r$metric), ]
  c$what <- sub(paste0("^", px, "_(pred|obs)_.*"), "\\1", c$metric); c$bin <- as.integer(sub(".*_", "", c$metric))
  w <- reshape(c[c("ID", "learner", "bin", "what", "value")], idvar = c("ID", "learner", "bin"), timevar = "what",
               direction = "wide")
  names(w) <- sub("^value\\.", "", names(w))
  q <- function(p) function(v) quantile(v, p, na.rm = TRUE)
  m <- aggregate(cbind(pred, obs) ~ learner + bin, w, mean)
  m$lo <- aggregate(obs ~ learner + bin, w, q(0.1))$obs; m$hi <- aggregate(obs ~ learner + bin, w, q(0.9))$obs
  m$model <- factor(lab(m$learner), lab(tb$learner)); m$family <- fam(m$learner); m
}
cs <- cal("calS")
f8 <- ggplot(cs, aes(pred, obs, colour = family, fill = family)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = 2) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.7) + geom_point(size = 1.3) +
  facet_wrap(~ model, nrow = 3) + scale_colour_manual(values = FAM_COL, name = NULL) +
  scale_fill_manual(values = FAM_COL, name = NULL) + coord_equal(xlim = c(0.4, 1), ylim = c(0.4, 1)) +
  labs(x = "Predicted event-free survival at 2 years (decile mean)", y = "Observed (Kaplan-Meier) in held-out third",
       title = "Survival calibration at 2 years, ACTG 175 (band: 10-90% of splits)") + th
save(f8, "fig8_actg175_calibration", 10, 7)

cb <- cal("calB")
f9 <- ggplot(cb, aes(pred * 100, obs * 100, colour = family, fill = family)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = 2) + geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo * 100, ymax = hi * 100), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.7) +
  geom_point(size = 1.3) + facet_wrap(~ model, nrow = 3) + scale_colour_manual(values = FAM_COL, name = NULL) +
  scale_fill_manual(values = FAM_COL, name = NULL) +
  labs(x = "Predicted benefit of ZDV + ddI at 2 years (quintile mean, points)",
       y = "Observed benefit (KM difference, points)",
       title = "Benefit calibration at 2 years, ACTG 175 (band: 10-90% of splits)") + th
save(f9, "fig9_actg175_benefit_calibration", 10, 7)

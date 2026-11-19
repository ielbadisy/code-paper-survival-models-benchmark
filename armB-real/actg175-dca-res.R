# armB-real/actg175-dca-res.R -- decision curves on ACTG 175 over 50 splits.
# In: results/raw_actg175_dca.rds (actg175-dca.R)
# Out: results/figures/fig10_actg175_dca (treatment), figS_actg175_nb (prognostic) .pdf/.png, dca_actg175.csv
# Run from armB-real/:  Rscript actg175-dca-res.R
suppressMessages(library(ggplot2))
R <- readRDS("results/raw_actg175_dca.rds"); r <- R$raw; OUT <- "results/figures"
COX <- c("coxph", "cox_int", "glmnet", "coxboost", "rpsurv", "rpsurv_tve")
LAB <- c(coxph = "Cox PH", cox_int = "Cox, all trt interactions", glmnet = "Ridge/lasso Cox", coxboost = "CoxBoost",
         rpsurv = "Royston-Parmar", rpsurv_tve = "Royston-Parmar TVE", nbsurv = "Naive Bayes",
         fastgbm = "GBM (fastgbm)", densemlp = "Neural net (MLP)", bart = "BART", xgboost = "XGBoost",
         rsf = "Random survival forest", orsf = "Oblique RSF", blackboost = "Tree boosting (mboost)")
FAM_COL <- c("Regression" = "#2166AC", "Machine learning" = "#D6604D")
fam <- function(l) factor(ifelse(l %in% COX, "Regression", "Machine learning"), names(FAM_COL))
th <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 11), plot.title.position = "plot")
save <- function(p, f, w, h) {
  ggsave(file.path(OUT, paste0(f, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT, paste0(f, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
}
q <- function(p) function(v) quantile(v, p)
summ <- function(px) {
  d <- r[grepl(px, r$metric), ]; d$x <- as.integer(sub(".*_", "", d$metric))
  m <- aggregate(value ~ learner + metric + x, d, mean)
  m$lo <- aggregate(value ~ learner + metric + x, d, q(0.1))$value
  m$hi <- aggregate(value ~ learner + metric + x, d, q(0.9))$value
  m$pos <- aggregate(value ~ learner + metric + x, d, function(v) mean(v > 0))$value; m
}

# ---- treatment decision curve: gain against delta, one panel per learner -------------------
g <- summ("^gain_"); write.csv(g, file.path(OUT, "dca_actg175.csv"), row.names = FALSE)
ord <- names(sort(tapply(g$value, g$learner, max), decreasing = TRUE))
g$model <- factor(LAB[g$learner], LAB[ord]); g$family <- fam(g$learner)
f10 <- ggplot(g, aes(x, 100 * value, colour = family, fill = family)) + geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_ribbon(aes(ymin = 100 * lo, ymax = 100 * hi), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.7) +
  facet_wrap(~ model, nrow = 3) + scale_colour_manual(values = FAM_COL, name = NULL) +
  scale_fill_manual(values = FAM_COL, name = NULL) +
  labs(x = "Cost of treatment (survival points at 2 years)", y = "Gain over treat-all / treat-none (points)") + th
save(f10, "fig10_actg175_dca", 10, 7)

# ---- prognostic decision curve: IPCW net benefit of flagging risk at 2 years ---------------
nb <- summ("^nb_"); na <- aggregate(value ~ x, summ("^nball_"), mean)
nb$family <- fam(nb$learner)
fS <- ggplot(nb, aes(x, value)) + geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(aes(group = learner, colour = family), linewidth = 0.5) +
  geom_line(data = na, linetype = 2, linewidth = 0.6) + scale_colour_manual(values = FAM_COL, name = NULL) +
  annotate("text", x = 22, y = -0.04, label = "Treat all", size = 3, hjust = 0) +
  coord_cartesian(ylim = c(-0.05, 0.25)) +
  labs(x = "Risk threshold at 2 years (%)", y = "Net benefit") + th
save(fS, "figS_actg175_nb", 6, 4.5)

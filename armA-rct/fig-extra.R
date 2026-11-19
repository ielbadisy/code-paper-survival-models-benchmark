# armA-rct/fig-extra.R -- whole-design figures for the Arm A v2 run: heatmap of differences from
# Cox PH over all scenario groups, nested loop plot of the 144 base cells (biplot is Fig 2A in fig-res.R).
# In: results/sim_raw_<tag>.rds. Out: results/figures_<tag>/fig{H,N}_*.{pdf,png}
# Run from armA-rct/:  Rscript fig-extra.R [tag]

suppressMessages({ library(ggplot2); library(patchwork) })
a <- commandArgs(TRUE); TAG <- if (length(a)) a[1] else "v2_90reps"
OUT <- sprintf("results/figures_%s", TAG)
COX <- c("coxph", "cox_spec", "glmnet", "coxboost", "rpsurv", "rpsurv_tve")
LAB <- c(coxph = "Cox PH", cox_spec = "Cox, true form", glmnet = "Ridge/lasso Cox", coxboost = "CoxBoost",
         rpsurv = "Royston-Parmar", rpsurv_tve = "Royston-Parmar TVE", nbsurv = "Naive Bayes",
         fastgbm = "GBM (fastgbm)", densemlp = "Neural net (MLP)", bart = "BART", xgboost = "XGBoost",
         rsf = "Random survival forest", orsf = "Oblique RSF", blackboost = "Tree boosting (mboost)")
FAM_COL <- c("Regression" = "#2166AC", "Machine learning" = "#D6604D")
HZ <- c(dec = "Decr.", const = "Const.", inc = "Incr.", nph = "Crossing")
ST <- c(nle = "Non-linear", tve = "Time-varying", censdep = "Dep. censoring")
fam <- function(l) factor(ifelse(l %in% COX, "Regression", "Machine learning"), names(FAM_COL))
lab <- function(l) unname(LAB[l])
th <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 11), plot.title.position = "plot")
save <- function(p, f, w, h) {
  ggsave(file.path(OUT, paste0(f, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT, paste0(f, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
}

r <- readRDS(sprintf("results/sim_raw_%s.rds", TAG))
cell <- c("dataset", "hazard", "beta_trt", "n", "cens", "stress")
M <- c("mae_true", "ibs", "cindex", "uno_c", "dwrange", "benefit_rmse", "benefit_cindex", "policy_regret")
r <- r[r$metric %in% M, ]
cm <- aggregate(value ~ ., r[c(cell, "learner", "metric", "value")], mean)
cm <- reshape(cm, idvar = c(cell, "learner"), timevar = "metric", direction = "wide")
names(cm) <- sub("^value\\.", "", names(cm))
# within-cell difference from Cox PH
ref <- cm[cm$learner == "coxph", c(cell, M)]; names(ref)[-(1:6)] <- paste0(M, "_ref")
cm <- merge(cm, ref, by = cell)
for (m in M) cm[[paste0(m, "_d")]] <- cm[[m]] - cm[[paste0(m, "_ref")]]

# ---- Heatmap: learner x scenario group, error and benefit RMSE vs Cox PH --------------------
dl <- ifelse(cm$dataset == "pbc", "PBC", "Prostate")
cm$block <- ifelse(cm$stress == "none", paste0(dl, "\n", HZ[cm$hazard]), "Stressors\n(n = 200)")
cm$col <- ifelse(cm$stress == "none", paste0("n=", cm$n), paste(dl, ST[cm$stress]))
hd <- aggregate(cbind(mae_true_d, benefit_rmse_d) ~ learner + block + col, cm[cm$learner != "coxph", ], mean)
hd$block <- factor(hd$block, c(paste0("PBC\n", HZ), paste0("Prostate\n", HZ), "Stressors\n(n = 200)"))
hd$col <- factor(hd$col, c(paste0("n=", c(100, 200, 400)), outer(c("PBC", "Prostate"), ST, paste)))
hd <- hd[order(hd$block, hd$col), ]
xm <- do.call(rbind, lapply(split(hd, hd$learner), function(d) c(d$mae_true_d, d$benefit_rmse_d)))
ord <- rownames(xm)[hclust(dist(xm))$order]
hd$model <- factor(lab(hd$learner), lab(ord))
L <- 0.04
hl <- function(v, tt) ggplot(hd, aes(col, model, fill = pmax(pmin(.data[[v]], L), -L))) +
  geom_tile(colour = "white", linewidth = 0.3) + facet_grid(~ block, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-L, L),
                       name = "Difference from Cox PH (blue = better, capped at \u00b10.04)  ",
                       guide = guide_colourbar(barwidth = unit(6, "cm"), title.vjust = 0.8)) +
  labs(x = NULL, y = NULL, title = tt) + th +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7), panel.grid = element_blank(),
        strip.text = element_text(size = 7.5), panel.spacing = unit(2, "pt"),
        axis.text.y = element_text(colour = FAM_COL[as.integer(fam(ord))]))
fH <- hl("mae_true_d", "A. Error against the true survival function") /
  hl("benefit_rmse_d", "B. Benefit RMSE") + plot_layout(guides = "collect") & theme(legend.position = "bottom")
save(fH, "figH_heatmap", 11, 7.5)

# ---- Nested loop plot: 144 base cells ------------------------------------------------------
b <- cm[cm$stress == "none", ]
FO <- c("dataset", "hazard", "beta_trt", "cens", "n")
b$hazard <- factor(b$hazard, names(HZ)); b$beta_trt <- factor(b$beta_trt, c(0, -0.4, 0.8))
key <- unique(b[FO]); key <- key[do.call(order, key), ]; key$x <- seq_len(nrow(key))
b <- merge(b, key, by = FO)
SHOW <- c("coxph", "glmnet", "cox_spec", "bart", "rsf", "xgboost")
SCOL <- setNames(c("black", "#2166AC", "#67A9CF", "#D6604D", "#B2182B", "#F4A582"), lab(SHOW))
b$model <- factor(lab(b$learner), c(lab(SHOW), setdiff(lab(names(LAB)), lab(SHOW))))
b$hi <- b$learner %in% SHOW
# step lines for the design factors, drawn under the curves
lv <- function(f) as.integer(factor(key[[f]], unique(sort(key[[f]]))))
nl <- function(v, yl, tt) {
  y0 <- min(b[[v]]); h <- diff(range(b[[v]])) * 0.05
  st <- do.call(rbind, lapply(seq_along(FO), function(i) {
    l <- lv(FO[i]); data.frame(x = key$x, y = y0 - h * (i * 2.2) + h * (l - 1) / max(l - 1, 1), f = FO[i])
  }))
  stl <- aggregate(y ~ f, st, min); stl$x <- -1
  ggplot(b, aes(x, .data[[v]], group = model)) +
    geom_step(data = b[!b$hi, ], colour = "grey80", linewidth = 0.3, direction = "mid") +
    geom_step(data = b[b$hi, ], aes(colour = model), linewidth = 0.55, direction = "mid") +
    geom_step(data = st, aes(x, y, group = f), inherit.aes = FALSE, linewidth = 0.4, direction = "mid") +
    geom_text(data = stl, aes(x, y, label = c(dataset = "Trial: PBC, prostate",
                                               hazard = "Hazard: decr., const., incr., crossing",
                                               beta_trt = "Log HR: 0, -0.4, 0.8", cens = "Censoring: 30%, 60%",
                                               n = "n: 100, 200, 400")[f]),
              inherit.aes = FALSE, hjust = 1, vjust = 0.2, size = 2.4) +
    scale_colour_manual(values = SCOL, name = NULL) + scale_x_continuous(expand = expansion(add = c(28, 1))) +
    scale_y_continuous(breaks = function(l) pretty(c(0, l[2]))) +
    labs(x = "Scenario (144 base scenarios, nested as in the steps below the curves)", y = yl, title = tt) + th +
    theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
}
fN <- nl("mae_true", "Mean error against the truth", "A. Survival prediction") /
  nl("benefit_rmse", "Benefit RMSE", "B. Treatment benefit") + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
save(fN, "figN_nested_loop", 11, 9)

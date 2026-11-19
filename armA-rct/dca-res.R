# armA-rct/dca-res.R -- treatment decision curves of the simulation (core scenarios, n = 200, 30% censoring):
# gain of "treat if predicted benefit > delta" over the better of treat-all and treat-none, against the
# optimal rule 1{b > delta} (headroom). Averaged over the four hazard shapes.
# In: results/dca_raw.rds (dca-run.R). Out: results/figures_v2_90reps/figD_dca.{pdf,png}, dca_data.csv
# Run from armA-rct/:  Rscript dca-res.R
suppressMessages(library(ggplot2))
OUT <- "results/figures_v2_90reps"
COX <- c("coxph", "cox_spec", "glmnet", "coxboost", "rpsurv", "rpsurv_tve")
LAB <- c(optimal = "Optimal rule (headroom)", coxph = "Cox PH", cox_spec = "Cox, true form", glmnet = "Ridge/lasso Cox",
         bart = "BART", rpsurv_tve = "Royston-Parmar TVE", rsf = "Random survival forest")
SCOL <- setNames(c("black", "grey30", "#67A9CF", "#2166AC", "#D6604D", "#053061", "#B2182B"), LAB)
th <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 11), plot.title.position = "plot")

r <- readRDS("results/dca_raw.rds"); r$delta <- as.integer(sub("gain_", "", r$metric))
by <- c("dataset", "hazard", "beta_trt", "learner", "delta")
m <- aggregate(value ~ ., r[c(by, "value")], mean)
m$se <- aggregate(value ~ ., r[c(by, "value")], function(v) sd(v) / sqrt(length(v)))$value
# average over hazards; MCSE of an average over C scenarios: sqrt(sum MCSE^2) / C
a <- aggregate(cbind(value, se2 = se^2) ~ dataset + beta_trt + learner + delta, m, mean)
a$mcse <- sqrt(a$se2 / 4); a$se2 <- NULL
write.csv(a, file.path(OUT, "dca_data.csv"), row.names = FALSE)
a$hi <- a$learner %in% names(LAB)
a$model <- factor(ifelse(a$hi, LAB[a$learner], "Other learners"), c(LAB, "Other learners"))
a$panel_r <- ifelse(a$dataset == "pbc", "PBC", "Prostate")
a$panel_c <- factor(paste("Log HR", a$beta_trt), paste("Log HR", c(0, -0.4, 0.8)))
h <- a[a$hi, ]
f <- ggplot(a, aes(delta, 100 * value, group = learner)) + geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(data = a[!a$hi, ], colour = "grey80", linewidth = 0.4) +
  geom_ribbon(data = h, aes(ymin = 100 * (value - 1.96 * mcse), ymax = 100 * (value + 1.96 * mcse), fill = model),
              alpha = 0.15, colour = NA) +
  geom_line(data = h, aes(colour = model, linetype = model), linewidth = 0.6) +
  facet_grid(panel_r ~ panel_c, scales = "free_y") +
  scale_colour_manual(values = SCOL, name = NULL) + scale_fill_manual(values = SCOL, guide = "none") +
  scale_linetype_manual(values = setNames(c(2, rep(1, 6)), LAB), name = NULL) +
  labs(x = expression("Cost of treatment " * delta * ", survival points at " * t^"*"),
       y = "Gain over treat-all / treat-none (points)") + th + guides(colour = guide_legend(nrow = 2))
ggsave(file.path(OUT, "figD_dca.pdf"), f, width = 9, height = 6)
ggsave(file.path(OUT, "figD_dca.png"), f, width = 9, height = 6, dpi = 200, bg = "white")

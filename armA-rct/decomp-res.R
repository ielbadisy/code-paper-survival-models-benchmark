# armA-rct/decomp-res.R -- summarise the benefit error decomposition (Appendix A, Proposition 1).
# In: results/analysis_v2_90reps/decomp.csv (from decomp-run.R).
# Out: results/analysis_v2_90reps/decomp_learner.csv, decomp_share.csv
# Run from armA-rct/:  Rscript decomp-res.R
P <- "results/analysis_v2_90reps/"
w <- read.csv(paste0(P, "decomp.csv"))
w$gap <- abs(w$bmse - 2 * w$pe * (1 - w$rho))   # identity check, should be ~0
lw <- aggregate(cbind(pe, rho, bmse, gap) ~ learner + dataset, w, mean)
write.csv(lw, paste0(P, "decomp_learner.csv"), row.names = FALSE)
# within each scenario, across learners: share of variance of log benefit MSE explained by log P alone,
# and Spearman of learner ranks by P and by benefit MSE
sh <- do.call(rbind, lapply(split(w, paste(w$dataset, w$hazard, w$beta_trt)), function(d)
  data.frame(d[1, c("dataset", "hazard", "beta_trt")], r2_pe = summary(lm(log(bmse) ~ log(pe), d))$r.squared,
             rho_rank = cor(d$pe, d$bmse, method = "spearman"), rho_sd = sd(d$rho))))
write.csv(sh, paste0(P, "decomp_share.csv"), row.names = FALSE)

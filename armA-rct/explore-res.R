# armA-rct/explore-res.R -- distil Arm A results: paired contrasts vs coxph,
# data-metric vs truth agreement (A3), prediction vs benefit, policy value, true benefit spread.
# In: results/sim_raw_<tag>.rds, results/design_<tag>.rds. Out: results/analysis_<tag>/*.csv
# Run from armA-rct/:  Rscript explore-res.R <tag>   (default v2_90reps)
source("../R/dgp.R"); source("../R/res-utils.R")
a <- commandArgs(TRUE); TAG <- if (length(a)) a[1] else "v2_90reps"
OUT <- sprintf("results/analysis_%s", TAG); dir.create(OUT, showWarnings = FALSE)
r <- readRDS(sprintf("results/sim_raw_%s.rds", TAG)); cell <- c("dataset","hazard","beta_trt","n","cens","stress")
W <- widen(r, c(cell, "ID", "learner"))
key <- do.call(paste, W[c(cell,"ID")])
out <- function(x, f) write.csv(x, file.path(OUT, f), row.names = FALSE)

# 1. paired difference vs coxph (same replicate), mean and MCSE, by dataset x hazard
M <- c("mae_true","ibs","dwrange","benefit_rmse","benefit_cindex","policy_regret","policy_gain")
cx <- W[W$learner == "coxph", ]; i <- match(key, do.call(paste, cx[c(cell,"ID")]))
D <- cbind(W[c(cell,"ID","learner")], W[M] - cx[i, M])
pd <- function(by) {
  a <- aggregate(D[M], D[by], function(v) mean(v, na.rm = TRUE))
  s <- aggregate(D[M], D[by], function(v) sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v))))
  names(s)[-seq_along(by)] <- paste0(M, "_mcse"); merge(a, s)
}
p1 <- pd(c("learner","dataset","hazard")); out(p1, "paired_vs_coxph_dataset_hazard.csv")
p2 <- pd(c("learner","n")); out(p2, "paired_vs_coxph_n.csv")
p3 <- pd(c("learner","stress")); out(p3, "paired_vs_coxph_stress.csv")

# 2. A3: within each cell-rep, does the learner an analyst picks on data metrics match the truth?
pick <- function(sc, best = min) do.call(rbind, lapply(split(W, key), function(d) {
  t <- d$mae_true; b <- d$benefit_rmse; s <- d[[sc]]; ok <- !is.na(s) & !is.na(t)
  j <- which(ok)[if (identical(best, min)) which.min(s[ok]) else which.max(s[ok])]
  data.frame(d[1, cell], picked = d$learner[j], mae_loss = t[j] - min(t, na.rm = TRUE),
             ben_loss = b[j] - min(b, na.rm = TRUE),
             rho = cor(if (identical(best, min)) s[ok] else -s[ok], t[ok], method = "spearman"))
}))
sel <- rbind(transform(pick("ibs"), by = "ibs"), transform(pick("dwrange"), by = "dwrange"),
             transform(pick("cindex", max), by = "harrell_c"), transform(pick("uno_c", max), by = "uno_c"))
out(sel, "selection_raw.csv")
out(aggregate(cbind(rho, mae_loss, ben_loss) ~ by + hazard, sel, mean), "selection_by_hazard.csv")
out(as.data.frame(prop.table(table(sel$by, sel$picked), 1)), "selection_picked.csv")

# 3. prediction vs benefit: per cell, Spearman between learner means of mae_true and benefit metrics
cm <- aggregate(W[M], W[c(cell,"learner")], mean, na.rm = TRUE)
pb <- do.call(rbind, lapply(split(cm, do.call(paste, cm[cell])), function(d)
  data.frame(d[1, cell], rho_rmse = cor(d$mae_true, d$benefit_rmse, method = "spearman"),
             rho_cben = cor(d$mae_true, -d$benefit_cindex, method = "spearman"),
             rho_regret = cor(d$mae_true, d$policy_regret, method = "spearman"))))
out(pb, "pred_vs_benefit_rank_cor.csv")

# 4. policy gain over the best of treat-all / treat-none, and benefit calibration slope (median)
pg <- aggregate(cbind(policy_gain, benefit_cal_slope, cal_slope_true) ~ learner + dataset + beta_trt + n, W,
                function(v) median(v), na.action = na.pass)
out(pg, "policy_gain_calslope.csv")

# 5. true benefit at t*: spread and share of patients harmed, per dataset x hazard x beta (stress none)
dz <- readRDS(sprintf("results/design_%s.rds", TAG))
g <- unique(dz$grid[dz$grid$stress == "none" & dz$grid$n == 200,
                    c("dataset", "hazard", "beta_trt", "cens", "cens_b")])
REFS <- list(pbc = load_reference("pbc"), prostate = load_reference("prostate"))
tb <- do.call(rbind, lapply(seq_len(nrow(g)), function(k) with(g[k, ], {
  ref <- REFS[[dataset]]; set.seed(1); X <- ref$X[sample.int(nrow(ref$X), 2e4, TRUE), , drop = FALSE]
  ts <- dz$times[[paste(dataset, hazard, beta_trt, "none", cens)]][5]
  b <- stress_true_surv(X, rep(1L, nrow(X)), ref, beta_trt, hazard, "none", ts) -
       stress_true_surv(X, rep(0L, nrow(X)), ref, beta_trt, hazard, "none", ts)
  data.frame(g[k, 1:4], sd_b = sd(b), q10 = quantile(b, .1), q90 = quantile(b, .9), harmed = mean(b < 0),
             headroom = min(mean(pmax(b, 0)), mean(pmax(-b, 0))))  # max gain of any rule over treat-all/none
})))
out(tb, "true_benefit_spread.csv")

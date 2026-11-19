# paper-numbers.R -- every result number quoted in the prose of the paper and of Supporting
# Information S1, recomputed from the stored outputs (no simulation is rerun). Run from the
# repository root: Rscript paper-numbers.R
# Inputs:  armA-rct/results/ (sim_raw_v2_90reps.rds, sim_summary_v2_90reps.csv, design_v2_90reps.rds,
#          tune_raw_v2.rds, dgp_check.csv, figures_v2_90reps/*.csv, analysis_v2_90reps/*.csv),
#          armB-real/results/ (raw_actg175_50splits.rds, figures/*.csv), data-ref/, BART::ACTG175
# Output:  results/paper_numbers.csv (quantity, value), one row per claim, in paper order.
# Scale: errors, MSEs and correlations as proportions; policy gains and headroom in survival points
# (x 100) where the text quotes points; shares as proportions.

suppressMessages(library(survival)); source("R/dgp.R")

A <- "armA-rct/results/"; FG <- paste0(A, "figures_v2_90reps/"); AN <- paste0(A, "analysis_v2_90reps/")
B <- "armB-real/results/"
ML <- c("nbsurv", "fastgbm", "densemlp", "bart", "xgboost", "rsf", "orsf", "blackboost")
REG5 <- c("coxph", "glmnet", "coxboost", "rpsurv", "rpsurv_tve")
SIGN <- c(mae_true = -1, ibs = -1, cindex = 1, uno_c = 1, dwrange = -1, benefit_rmse = -1,
          benefit_cindex = 1, policy_regret = -1)

out <- list()
put <- function(...) out <<- c(out, list(...))
rd <- function(f) read.csv(f, stringsAsFactors = FALSE)
at <- function(d, l, col) d[[col]][match(l, d$learner)]
rg <- function(k, x) setNames(range(x), paste0(k, c("_min", "_max")))
pts <- function(x) 100 * x

## Methods: design ---------------------------------------------------------------------------
ref <- lapply(c(pbc = "pbc", prostate = "prostate"), load_reference)
g <- readRDS(paste0(A, "design_v2_90reps.rds"))$grid
base <- g$stress == "none"; core <- base & g$n == 200 & g$cens == 0.3
put(pbc_n = nrow(ref$pbc$X), pbc_p = ncol(ref$pbc$X), pro_n = nrow(ref$prostate$X),
    pro_p = ncol(ref$prostate$X), n_base = sum(base), n_stress = sum(!base), n_scen = nrow(g),
    n_core = sum(core), n_matched = nrow(unique(g[!base, c("dataset", "hazard")])))

dg <- rd(paste0(A, "dgp_check.csv"))
put(rg("dgp_cens", dg$cens), dgp_truth_max = max(dg$truth_vs_latent),
    rg("dgp_km_censdep", dg$km_vs_truth[dg$stress == "censdep"]))

tb <- rd(paste0(AN, "true_benefit_spread.csv")); tp <- tb[tb$dataset == "prostate", ]
put(rg("pro_sd_b", tp$sd_b), pro_harmed_max_null = max(tp$harmed[tp$beta_trt == 0]),
    pro_harmed_max = max(tp$harmed))

tu <- readRDS(paste0(A, "tune_raw_v2.rds")); ng <- lengths(tu$hp)
put(n_tune_groups = nrow(unique(tu$best[c("dataset", "hazard", "n")])), rg("grid_size", ng[ng > 0]),
    n_learners = length(readRDS(paste0(A, "design_v2_90reps.rds"))$learners))

## Methods: MCSE, failures, software ---------------------------------------------------------
s <- rd(paste0(A, "sim_summary_v2_90reps.csv")); m <- s$mcse[s$metric == "mae_true"]
r <- readRDS(paste0(A, "sim_raw_v2_90reps.rds"))
ib <- r[r$metric == "ibs", ]; fl <- ib[is.na(ib$value), ]
se <- r[r$metric == "secs", ]; sl <- tapply(se$value, se$learner, mean)
put(mcse_mae_max = max(m), mcse_mae_median = median(m), n_fail = nrow(fl), n_fits = nrow(ib),
    fail_n = fl$n, fail_cens = pts(fl$cens), fail_beta = fl$beta_trt,
    secs_mean = mean(se$value), secs_coxph = sl[["coxph"]], secs_bart = sl[["bart"]])
v <- readRDS(paste0(A, "design_v2_90reps.rds"))$versions
put(as.list(setNames(as.character(v), paste0("ver_", names(v)))))

## Methods: ACTG 175 -------------------------------------------------------------------------
data(ACTG175, package = "BART")
a <- subset(ACTG175, arms %in% 0:1); a$trt <- as.integer(a$arms == 1); a$time <- a$days / 365.25
km <- summary(survfit(Surv(time, cens) ~ trt, a), times = 2)$surv
# "prior antiretroviral therapy" in the text is str2 (antiretroviral history)
pint <- summary(coxph(Surv(time, cens) ~ trt * str2, a))$coef["trt:str2", 5]
put(actg_n0 = sum(a$trt == 0), actg_n1 = sum(a$trt == 1), actg_avg_ben = pts(diff(km)))

## Results: survival prediction --------------------------------------------------------------
t1 <- rd(paste0(FG, "tab1_core.csv")); tf <- rd(paste0(FG, "tab1_full.csv"))
for (l in c("cox_spec", "rpsurv_tve", "coxboost", "glmnet", "coxph", "rpsurv", "bart", "rsf", "densemlp",
            "xgboost", "nbsurv")) out[[paste0("core_mae_", l)]] <- at(t1, l, "mae_true")
put(rg("core_mae_mcse", t1$mae_true_mcse))
rc <- r[r$metric == "mae_true" & r$stress == "none" & r$n == 200 & r$cens == 0.3, ]
mh <- tapply(rc$value, list(rc$learner, rc$hazard), mean); dh <- sweep(mh, 2, mh["coxph", ])
put(rg("nph_gain_tve", -dh[c("cox_spec", "rpsurv_tve"), "nph"]),
    rg("other_gain_shrink", -dh[c("coxboost", "glmnet"), c("dec", "const", "inc")]),
    full_mae_coxboost = at(tf, "coxboost", "mae_true"), full_rank_coxboost = at(tf, "coxboost", "rank_pred"),
    full_mae_glmnet = at(tf, "glmnet", "mae_true"), full_rank_glmnet = at(tf, "glmnet", "rank_pred"),
    full_mae_oracle = at(tf, "cox_spec", "mae_true"), full_rank_oracle = at(tf, "cox_spec", "rank_pred"))

f1 <- rd(paste0(FG, "fig1_data.csv")); n1 <- f1[f1$n == 100, ]; n4 <- f1[f1$n == 400, ]
put(n100_gain_coxboost = -at(n1, "coxboost", "mae_true"), n100_gain_glmnet = -at(n1, "glmnet", "mae_true"),
    rg("n100_gain_orsf_bart", -at(n1, c("orsf", "bart"), "mae_true")),
    n400_shrink_absdiff_max = max(abs(at(n4, c("glmnet", "coxboost"), "mae_true"))),
    n400_ml_worse = sum(at(n4, ML, "mae_true") > 0), n400_loss_bart = at(n4, "bart", "mae_true"),
    n400_loss_orsf = at(n4, "orsf", "mae_true"),
    rg("n400_loss_mlp_rsf", at(n4, c("densemlp", "rsf"), "mae_true")),
    n400_loss_xgboost = at(n4, "xgboost", "mae_true"),
    rsf_worse_all_n = sum(f1$mae_true[f1$learner == "rsf"] > 0))

## Results: choosing a learner from the data -------------------------------------------------
t2 <- rd(paste0(FG, "tab2_core.csv")); b2 <- function(k, col) t2[[col]][t2$by == k]
put(sel_rho_ibs = b2("ibs", "rho"), sel_loss_ibs = b2("ibs", "mae_loss"),
    sel_rho_harrell = b2("harrell_c", "rho"), sel_loss_harrell = b2("harrell_c", "mae_loss"),
    sel_rho_harrell_nph = b2("harrell_c", "rho_nph"), sel_loss_harrell_nph = b2("harrell_c", "loss_nph"),
    sel_ratio_nph = b2("harrell_c", "loss_nph") / b2("ibs", "loss_nph"),
    sel_rho_uno = b2("uno_c", "rho"), sel_rho_dw = b2("dwrange", "rho"),
    sel_rho_mcse_max = max(t2$rho_mcse), sel_loss_mcse_max = max(t2$loss_mcse))

## Results: individual treatment benefit -----------------------------------------------------
bv <- rd(paste0(FG, "biplot_variance.csv"))
rk <- r[r$metric %in% names(SIGN) & r$stress == "none" & r$n == 200 & r$cens == 0.3, ]
bm <- tapply(rk$value, list(rk$learner, rk$metric), mean, na.rm = TRUE)[, names(SIGN)]
cm <- cor(sweep(bm, 2, SIGN, `*`)); pr <- names(SIGN)[1:5]
put(biplot_var2 = sum(bv$var_explained[1:2]), rg("cor_pred", cm[pr, pr][upper.tri(cm[pr, pr])]),
    rg("cor_ben_pred", cm[c("benefit_rmse", "benefit_cindex"), pr]),
    rg("cor_regret_pred", cm["policy_regret", pr]),
    rank_pred_coxboost = at(t1, "coxboost", "rank_pred"), rank_pred_glmnet = at(t1, "glmnet", "rank_pred"),
    rank_ben_coxboost = at(t1, "coxboost", "rank_benefit"),
    rank_ben_glmnet = at(t1, "glmnet", "rank_benefit"))
for (l in c("bart", "rpsurv_tve", "cox_spec", "coxph", "rsf"))
  out[[paste0("core_brmse_", l)]] <- at(t1, l, "benefit_rmse")
mid <- t1$benefit_rmse[t1$rank_benefit %in% 4:11]
put(n_brmse_mid = length(mid), rg("brmse_mid", mid), brmse_mcse_max = max(t1$benefit_rmse_mcse),
    cben_oracle = at(t1, "cox_spec", "benefit_cindex"), cben_bart = at(t1, "bart", "benefit_cindex"),
    cben_coxph = at(t1, "coxph", "benefit_cindex"))

dl <- rd(paste0(AN, "decomp_learner.csv")); ds <- rd(paste0(AN, "decomp_share.csv"))
dp <- dl[dl$dataset == "prostate", ]; db <- dl[dl$dataset == "pbc", ]; n0 <- ds[ds$beta_trt == 0, ]
put(decomp_r2_mean = mean(ds$r2_pe), rg("decomp_r2_null", tapply(n0$r2_pe, n0$dataset, mean)),
    rg("pro_rho_shrink", at(dp, c("glmnet", "coxboost"), "rho")), pro_rho_oracle = at(dp, "cox_spec", "rho"),
    pro_bmse_ratio = at(dp, "cox_spec", "bmse") / at(dp, "glmnet", "bmse"))

f2 <- rd(paste0(FG, "fig2b_medians.csv")); fn <- f2[f2$stress == "none", ]
put(agree_neg = fn$rho_rmse[fn$beta_trt == -0.4], agree_pos = fn$rho_rmse[fn$beta_trt == 0.8],
    agree_null = fn$rho_rmse[fn$beta_trt == 0], rg("agree_stress", f2$rho_rmse[f2$stress != "none"]))
for (l in c("coxph", "bart", "fastgbm", "orsf", "rsf"))
  out[[paste0("bslope_", l)]] <- at(t1, l, "benefit_cal_slope")

f4 <- rd(paste0(FG, "fig4_data.csv")); st <- function(k) f4[f4$stress == k, ]
T3 <- c("cox_spec", "bart", "rpsurv_tve")
put(stress_none_oracle = -at(st("none"), "cox_spec", "benefit_rmse"),
    stress_none_bart = -at(st("none"), "bart", "benefit_rmse"),
    stress_none_rptve = -at(st("none"), "rpsurv_tve", "benefit_rmse"),
    stress_none_top3 = sum(rank(st("none")$benefit_rmse)[match(T3, st("none")$learner)] <= 3),
    stress_tve_oracle = -at(st("tve"), "cox_spec", "benefit_rmse"),
    stress_tve_bart = -at(st("tve"), "bart", "benefit_rmse"),
    stress_tve_rptve = -at(st("tve"), "rpsurv_tve", "benefit_rmse"),
    rg("stress_nle_top3", -at(st("nle"), T3, "benefit_rmse")),
    rg("stress_censdep_top3", -at(st("censdep"), T3, "benefit_rmse")),
    rg("stress_rsf_loss", f4$benefit_rmse[f4$learner == "rsf" & f4$stress != "censdep"]))

## Results: treatment rules ------------------------------------------------------------------
pg <- r[r$metric == "policy_gain" & r$dataset == "pbc", ]
f5 <- rd(paste0(FG, "fig5_data.csv")); be <- function(b) f5[f5$beta_trt == b, ]
hd <- tapply(f5$headroom, f5$beta_trt, mean)
put(pbc_gain_max = max(pg$value, na.rm = TRUE), headroom_null = pts(hd[["0"]]),
    headroom_neg = pts(hd[["-0.4"]]), headroom_pos = pts(hd[["0.8"]]),
    rg("neg_loss", -pts(be(-0.4)$policy_gain)), null_gain_oracle = pts(at(be(0), "cox_spec", "policy_gain")),
    rg("null_gain_mlp_bart_rsf", pts(at(be(0), c("densemlp", "bart", "rsf"), "policy_gain"))),
    null_loss_coxph = -pts(at(be(0), "coxph", "policy_gain")),
    pos_n_gain_01 = sum(be(0.8)$policy_gain > 0.001),
    pos_gain_oracle = pts(at(be(0.8), "cox_spec", "policy_gain")),
    regret_oracle = at(t1, "cox_spec", "policy_regret"), regret_bart = at(t1, "bart", "policy_regret"),
    regret_max = max(t1$policy_regret),
    regret_max_is_blackboost = as.numeric(t1$learner[which.max(t1$policy_regret)] == "blackboost"))

dc <- rd(paste0(FG, "dca_data.csv")); dc <- dc[dc$beta_trt == -0.4, ]
dd <- function(ds, l) dc[dc$dataset == ds & dc$learner == l, ]
# PBC learners are read at the cost where the optimal rule peaks
op <- dd("pbc", "optimal"); dx <- op$delta[which.max(op$value)]
at8 <- function(l) with(dd("pbc", l), pts(value[delta == dx]))
oq <- dd("prostate", "optimal")
oth <- dc[dc$dataset == "prostate" & !dc$learner %in% c("optimal", "cox_spec", "bart"), ]
put(dca_pbc_opt0 = pts(op$value[op$delta == 0]), dca_pbc_opt_max = pts(max(op$value)), dca_pbc_opt_at = dx,
    dca_pbc_rptve = at8("rpsurv_tve"), dca_pbc_oracle = at8("cox_spec"),
    rg("dca_pbc_loss_coxph_glmnet", -c(at8("coxph"), at8("glmnet"))),
    dca_pro_opt_max = pts(max(oq$value)), dca_pro_opt_at = oq$delta[which.max(oq$value)],
    dca_pro_oracle_max = pts(max(dd("prostate", "cox_spec")$value)),
    dca_pro_bart_max = pts(max(dd("prostate", "bart")$value)), dca_pro_other_max = pts(max(oth$value)),
    dca_pro_rsf_min = pts(min(dd("prostate", "rsf")$value)),
    dca_mcse_max = pts(max(dc$mcse)))

## Results: ACTG 175 -------------------------------------------------------------------------
x <- readRDS(paste0(B, "raw_actg175_50splits.rds"))$raw; ta <- rd(paste0(B, "figures/tab_actg175.csv"))
nl <- rd(paste0(B, "figures/null_actg175.csv")); WORSE <- c("densemlp", "nbsurv", "cox_int")
mt <- function(k) tapply(x$value[x$metric == k], x$learner[x$metric == k], mean)
ii <- mt("ibs"); ic <- mt("cindex"); good <- setdiff(names(ii), WORSE)
isd <- tapply(x$value[x$metric == "ibs"], x$learner[x$metric == "ibs"], sd)
put(actg_km_ctrl = pts(km[1]), actg_km_trt = pts(km[2]),
    actg_hr = exp(coef(coxph(Surv(time, cens) ~ trt, a)))[[1]], actg_fail = sum(ta$fail),
    actg_n_close = length(good), rg("actg_ibs_close", ii[good]), actg_ibs_km = nl$ibs[1],
    actg_ibs_km_arm = nl$ibs[2],
    actg_ibs_gain_max = nl$ibs[1] - min(ii), rg("actg_ibs_sd", isd[good]), rg("actg_ibs_worse", ii[WORSE]),
    actg_rho_ibs_harrell = cor(rank(ii), rank(-ic), method = "spearman"),
    actg_rank_ibs_coxph = rank(ii)[["coxph"]], actg_rank_harrell_coxph = rank(-ic)[["coxph"]])

BEN <- c("coxph", "rpsurv", "glmnet", "coxboost", "blackboost", "fastgbm", "xgboost")
q <- x[x$metric %in% c("calB_pred_1", "calB_pred_5") & x$learner %in% BEN, ]
qm <- tapply(q$value, list(q$learner, q$metric), mean)
sp <- x[x$metric == "benefit_spread", ]; sp <- sp[order(sp$ID), ]; sx <- function(l) sp$value[sp$learner == l]
cx <- ta[ta$learner == "cox_int", ]
put(rg("actg_pred_ben", pts(c(qm[, 1], qm[, 2]))), rg("actg_bslope_reg", at(ta, BEN, "benefit_cal_slope")),
    actg_spread_rsf = pts(at(ta, "rsf", "benefit_spread")),
    actg_spread_orsf = pts(at(ta, "orsf", "benefit_spread")),
    actg_spread_coxph = pts(at(ta, "coxph", "benefit_spread")),
    rg("actg_forest_beat", c(mean(sx("rsf") > sx("coxph")), mean(sx("orsf") > sx("coxph")))),
    actg_quintile_arm = mean(table(a$trt)) / 3 / 5, actg_p_int = pint,
    actg_coxint_withheld = pts(1 - cx$treat_share), actg_coxint_loss = -pts(cx$policy_gain))

gd <- rd(paste0(B, "figures/dca_actg175.csv")); w <- tapply(gd$value, list(gd$learner, gd$x), mean)
p8 <- gd$pos[gd$x == 8 & gd$learner %in% c("glmnet", "coxboost")]
put(rg("actg_dca_reg_6_8", pts(w[REG5, as.character(6:8)])), actg_dca_pos8 = min(p8),
    actg_dca_ml_max = pts(max(w[ML, ])),
    actg_dca_n_nonlosing_above12 = max(colSums(w[, as.numeric(colnames(w)) > 12] >= 0)),
    sim_brmse_rsf = at(t1, "rsf", "benefit_rmse"), sim_brmse_orsf = at(t1, "orsf", "benefit_rmse"),
    sim_cben_rsf = at(t1, "rsf", "benefit_cindex"), sim_cben_orsf = at(t1, "orsf", "benefit_cindex"))

## Discussion --------------------------------------------------------------------------------
ps <- rd(paste0(AN, "paired_vs_coxph_stress.csv")); ps <- ps[ps$stress != "none", ]
pm <- function(l, k) ps$mae_true[ps$learner == l & ps$stress %in% k]
put(rg("stress_gain_tve2", -c(pm("cox_spec", ps$stress), pm("rpsurv_tve", ps$stress))),
    stress_gain_bart_nle = -pm("bart", "nle"),
    stress_gain_shrink_max = max(-c(pm("glmnet", ps$stress), pm("coxboost", ps$stress))),
    stress_loss_glmnet_nle = pm("glmnet", "nle"))

## Appendix A --------------------------------------------------------------------------------
hb <- function(b) tp$headroom[tp$beta_trt == b]
put(rg("H_neg", hb(-0.4)), H_pos_max = max(hb(0.8)), rg("H_null", hb(0)),
    rg("pro_pe_shrink", at(dp, c("glmnet", "coxboost"), "pe")),
    pro_bmse_shrink = max(at(dp, c("glmnet", "coxboost"), "bmse")),
    pro_pe_oracle = at(dp, "cox_spec", "pe"), pro_bmse_oracle = at(dp, "cox_spec", "bmse"),
    pbc_pe_bart = at(db, "bart", "pe"), pbc_pe_coxph = at(db, "coxph", "pe"),
    pbc_rho_bart = at(db, "bart", "rho"), pbc_rho_coxph = at(db, "coxph", "rho"),
    pbc_bmse_bart = at(db, "bart", "bmse"), pbc_bmse_coxph = at(db, "coxph", "bmse"),
    pbc_null_rank_rho = mean(n0$rho_rank[n0$dataset == "pbc"]), decomp_gap_max = max(abs(dl$gap)))

o <- unlist(out)
dir.create("results", showWarnings = FALSE)
write.csv(data.frame(quantity = names(o), value = unname(o)), "results/paper_numbers.csv", row.names = FALSE)

# R/res-utils.R -- shared by the *-res.R scripts: long-to-wide results, learner labels and families,
# plot theme, figure saving (pdf + png into OUT, set by the calling script).
suppressMessages({ library(ggplot2); library(patchwork) })
COX <- c("coxph", "cox_spec", "glmnet", "coxboost", "rpsurv", "rpsurv_tve", "cox_int")
LAB <- c(coxph = "Cox PH", cox_spec = "Cox, true form", glmnet = "Ridge/lasso Cox", coxboost = "CoxBoost",
         rpsurv = "Royston-Parmar", rpsurv_tve = "Royston-Parmar TVE", nbsurv = "Naive Bayes",
         fastgbm = "GBM (fastgbm)", densemlp = "Neural net (MLP)", bart = "BART", xgboost = "XGBoost",
         rsf = "Random survival forest", orsf = "Oblique RSF", blackboost = "Tree boosting (mboost)",
         cox_int = "Cox, all trt interactions")
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
# one column per metric, one row per id combination
widen <- function(d, id, tv = "metric") {
  w <- reshape(d[c(id, tv, "value")], idvar = id, timevar = tv, direction = "wide")
  names(w) <- sub("^value\\.", "", names(w)); w
}

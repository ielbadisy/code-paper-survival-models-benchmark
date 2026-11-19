# armB-real/actg175-null.R -- covariate-free reference for the ACTG 175 IBS: Kaplan-Meier pooled and by arm,
# scored like the learners (same times, same IBS code) on 50 random 2:1 splits.
# In: results/raw_actg175_50splits.rds (times). Out: results/figures/null_actg175.csv
# Run from armB-real/:  Rscript actg175-null.R
setwd("../armA-rct"); suppressMessages(source("common.R")); setwd("../armB-real")
data(ACTG175, package = "BART"); A <- ACTG175[ACTG175$arms %in% 0:1, ]
D <- data.frame(time = A$days / 365.25, status = A$cens, trt = as.integer(A$arms == 1))
TIMES <- readRDS("results/raw_actg175_50splits.rds")$times
kmt <- function(d) summary(survfit(Surv(time, status) ~ 1, d), times = TIMES, extend = TRUE)$surv
ibs <- function(S, y) evaluate(S, y = y, times = TIMES, metrics = "ibs")$summary$mean
set.seed(2028)
v <- t(replicate(50, {
  i <- unlist(lapply(split(seq_len(nrow(D)), D$trt), function(k) sample(k, round(length(k) / 3))))
  tr <- D[-i, ]; te <- D[i, ]; y <- Surv(te$time, te$status); ka <- sapply(0:1, function(a) kmt(tr[tr$trt == a, ]))
  c(km = ibs(matrix(kmt(tr), nrow(te), length(TIMES), byrow = TRUE), y), km_arm = ibs(t(ka[, te$trt + 1]), y))
}))
write.csv(data.frame(model = c("Kaplan-Meier", "Kaplan-Meier by arm"), ibs = colMeans(v), sd = apply(v, 2, sd)),
          "results/figures/null_actg175.csv", row.names = FALSE)

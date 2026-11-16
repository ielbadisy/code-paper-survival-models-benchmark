# armA-rct/sim-res.R -- summarise Arm A raw results: mean and MCSE per
# cell x learner x metric, per-cell ranks, mean ranks, failure rates.
# Run from armA-rct/:  Rscript sim-res.R <tag>   (tag: smoke | 200reps ...)

args <- commandArgs(TRUE)
tag <- if (length(args)) args[1] else "200reps"
raw <- read.csv(sprintf("results/sim_raw_%s.csv", tag), stringsAsFactors = FALSE)
cell <- c("dataset", "hazard", "beta_trt", "n", "cens", "stress")

# lower is better for these, higher for the rest
lower_better <- c("ibs", "iae", "ise", "brier", "ece", "mae_true", "mae_true_tstar",
                  "benefit_rmse", "policy_regret", "secs")

agg <- aggregate(value ~ ., data = raw[c(cell, "learner", "metric", "value")],
                 FUN = function(v) c(mean = mean(v), median = stats::median(v),  # median: benefit_cal_slope is heavy-tailed
                                     mcse = stats::sd(v) / sqrt(length(v)), reps = length(v)))
summ <- cbind(agg[setdiff(names(agg), "value")], as.data.frame(agg$value))
fails <- aggregate(value ~ ., data = transform(raw[raw$metric == "ibs", c(cell, "learner", "value")],
                                               value = is.na(value)), FUN = mean, na.action = na.pass)
names(fails)[names(fails) == "value"] <- "fail_rate"

# per-cell rank of each learner on each metric (1 = best)
summ$rank <- ave(ifelse(summ$metric %in% lower_better, summ$mean, -summ$mean),
                 interaction(summ[c(cell, "metric")], drop = TRUE),
                 FUN = function(x) rank(x, na.last = "keep", ties.method = "average"))
mean_rank <- aggregate(rank ~ learner + metric, data = summ, FUN = mean)
mean_rank <- mean_rank[order(mean_rank$metric, mean_rank$rank), ]

dir.create("results", showWarnings = FALSE)
write.csv(summ, sprintf("results/sim_summary_%s.csv", tag), row.names = FALSE)
write.csv(fails, sprintf("results/sim_fail_rates_%s.csv", tag), row.names = FALSE)
write.csv(mean_rank, sprintf("results/sim_mean_rank_%s.csv", tag), row.names = FALSE)

# timing per learner (seconds per fit + 3 predictions), averaged over cells
timing <- aggregate(mean ~ learner + n, data = summ[summ$metric == "secs", ], FUN = mean)
write.csv(timing, sprintf("results/timing_%s.csv", tag), row.names = FALSE)
print(timing[order(-timing$mean), ], row.names = FALSE)

# armA-rct/tune-run.R -- learner tuning for the v2 benchmark (replicates 1-10).
# Per dataset x hazard x n (24 groups), 10 generated training trials (beta
# -0.4, 30% censoring, no stress). Each config in HP is fit on the training
# trial and scored by IBS on an independent validation trial (N = 500) from
# the same DGP; the config with the best mean IBS over the 10 is frozen.
# Output: results/tuned_v2.rds  TUNED[[dataset]][[hazard]][[n]][[learner]]
#         results/tune_raw_v2.rds (every config x replicate)
# Run from armA-rct/:  OMP_NUM_THREADS=1 Rscript tune-run.R <cores> [reps]

source("common.R")

xg <- function(...) { g <- expand.grid(..., stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  lapply(seq_len(nrow(g)), function(i) as.list(g[i, , drop = FALSE])) }
HP <- list(
  glmnet = xg(alpha = c(0, 0.5, 1), s = c("lambda.min", "lambda.1se")),
  coxboost = xg(stepno = c(25, 50, 100, 200), pmult = c(4.5, 9)),
  nbsurv = xg(cov_structure = c("diagonal", "full"), time_smooth = c(FALSE, TRUE)),
  fastgbm = xg(ntrees = c(50, 100, 200, 400), learning_rate = 0.03, max_depth = c(1, 2)),
  densemlp = with(expand.grid(h = 1:2, nb = c(5, 10), dr = c(0.1, 0.3)),
    Map(function(h, nb, dr) list(hidden_units = list(16, c(32, 16))[[h]], n_bins = nb, dropout = dr), h, nb, dr)),
  bart = xg(K = c(5, 10), power = c(1, 2), base = c(0.75, 0.95)),
  xgboost = xg(nrounds = c(50, 150), max_depth = c(2, 4), eta = c(0.05, 0.2)),
  rsf = xg(mtry = c(2, 4), nodesize = c(5, 15, 30)),
  orsf = xg(mtry = c(2, 4), leaf_min_obs = c(5, 15, 30)),
  blackboost = xg(mstop = c(50, 100, 200), maxdepth = c(1, 2)))

ibs_of <- function(S, y, times) evaluate(S, y = y, times = times, metrics = "ibs")$summary$mean

tune_one <- function(dataset, hazard, n, beta_trt, cens, stress, cens_b) {
  ref <- REFS[[dataset]]; times <- TIMES[[cell_key(dataset, hazard, beta_trt, stress, cens)]]
  train <- simulate_stress(ref, n, beta_trt, hazard, stress, cens_b); attr(train, "truth") <- NULL
  val <- simulate_stress(ref, N_TEST, beta_trt, hazard, stress, cens_b); attr(val, "truth") <- NULL
  y <- Surv(val$time, val$status)
  do.call(rbind, lapply(names(HP), function(l) do.call(rbind, lapply(seq_along(HP[[l]]), function(k)
    data.frame(learner = l, config = k, ibs = tryCatch({
      S <- predict_any(fit_any(l, train, ref, HP[[l]][[k]]), val, times)
      if (anyNA(S)) NA_real_ else ibs_of(S, y, times)
    }, error = function(e) NA_real_))))))
}

if (sys.nframe() == 0L) {
  a <- commandArgs(TRUE)
  cores <- if (length(a) >= 1) as.integer(a[1]) else 2L
  reps <- if (length(a) >= 2) as.integer(a[2]) else 10L
  tg <- unique(GRID[GRID$stress == "none", c("dataset", "hazard", "n")])
  tg <- prepare_grid(transform(tg, beta_trt = -0.4, cens = 0.3, stress = "none"))
  rownames(tg) <- NULL
  TIMES <- grid_times(tg)
  tag <- if (reps == 10L) "v2" else sprintf("v2_%dreps", reps)
  raw <- runsim(n = reps, grid_params = tg, sim_func = tune_one, seed = 1027, num_cores = cores,
                show_progress = FALSE, backend = "fork", on_error = "warn",
                save_path = sprintf("results/tune_checkpoint_%s.rds", tag), checkpoint_every = 1L)
  agg <- aggregate(ibs ~ dataset + hazard + n + learner + config, raw, mean)
  best <- do.call(rbind, lapply(split(agg, agg[c("dataset", "hazard", "n", "learner")], drop = TRUE),
                                function(d) d[which.min(d$ibs), ]))
  tuned <- list()
  for (i in seq_len(nrow(best))) with(best[i, ],
    tuned[[dataset]][[hazard]][[as.character(n)]][[learner]] <<- HP[[learner]][[config]])
  saveRDS(list(raw = raw, best = best, hp = HP, versions = pkg_versions()), sprintf("results/tune_raw_%s.rds", tag))
  saveRDS(tuned, sprintf("results/tuned_%s.rds", tag))
}

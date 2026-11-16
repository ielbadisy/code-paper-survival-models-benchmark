# armA-rct/tune-run.R -- tune each tunable Arm A learner ONCE per dataset x
# hazard (PROTOCOL.md Section 10, decided 2026-11-15), then freeze the
# settings for every replicate and cell of that dataset x hazard.
# Tuning trial: N = 1000, beta_trt = -0.4, 30% censoring, no stressor.
# Budget: survalis tune(), 3-fold CV, IBS, at most MAX_CONFIGS configurations
# of the learner's built-in grid. Same budget for every learner.
# glmnet is tuned with trt unpenalised (as in sim-run.R); coxboost is tuned
# through survalis (trt penalised) and the chosen stepno/penalty are reused
# by the unpenalised-trt fit in sim-run.R.
# Untuned (no built-in grid): coxph, cox_spec, rpsurv, rpsurv_tve, fastgbm, nbsurv.
# flexsurvreg stays at its Weibull default: its built-in tuner errors when one
# candidate distribution fails to converge, and the distribution is a model
# choice rather than a hyperparameter.
# Run from armA-rct/:  Rscript tune-run.R <cores>

suppressMessages({ library(survalis); library(survival) })
source("../R/dgp-rct-stress.R")

MAX_CONFIGS <- 8L
N_TUNE <- 1000L
TUNABLE <- c("glmnet", "coxboost", "blackboost", "rsf", "ranger",
             "orsf", "cforest", "survsvm")
metric_ids <- listmetrics()$metric

capped_grid <- function(learner) {
  f <- get(paste0("tune_", learner), envir = asNamespace("survalis"))
  g <- tryCatch(eval(formals(f)$param_grid, envir = environment(f)), error = function(e) NULL)
  if (!is.data.frame(g) || nrow(g) <= MAX_CONFIGS) return(NULL)  # built-in grid as is
  g[unique(round(seq(1, nrow(g), length.out = MAX_CONFIGS))), , drop = FALSE]
}

tune_task <- function(dataset, hazard, learner) {
  ref <- load_reference(dataset)
  set.seed(2026)
  b <- calibrate_stress_censoring(ref, -0.4, hazard, "none", 0.3, n_pilot = 5e4)
  d <- simulate_stress(ref, N_TUNE, -0.4, hazard, "none", b); attr(d, "truth") <- NULL
  times <- unname(quantile(d$time[d$status == 1], seq(0.1, 0.9, 0.1)))
  args <- list(Surv(time, status) ~ ., d, model = learner, grid = capped_grid(learner),
               times = times, resampling = cv(v = 3, seed = 1), metric = "ibs")
  if (learner == "glmnet") args$penalty.factor <- as.numeric(setdiff(names(d), c("time", "status")) != "trt")
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(do.call(tune, args), error = function(e) conditionMessage(e))
  secs <- proc.time()[["elapsed"]] - t0
  if (is.character(res)) return(list(dataset = dataset, hazard = hazard, learner = learner,
                                     params = NULL, secs = secs, error = res))
  best <- res$best
  params <- as.list(best[setdiff(names(best), c(metric_ids, "failed"))])
  list(dataset = dataset, hazard = hazard, learner = learner, params = params,
       ibs = best$ibs, secs = secs, error = NA)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  cores <- if (length(args)) as.integer(args[1]) else 2L
  tasks <- expand.grid(learner = TUNABLE, hazard = c("dec", "const", "inc", "nph"),
                       dataset = c("pbc", "prostate"), stringsAsFactors = FALSE)
  # slow learners first so the pool drains evenly
  tasks <- tasks[order(!tasks$learner %in% c("cforest", "survsvm", "rsf", "ranger")), ]
  out <- parallel::mclapply(seq_len(nrow(tasks)), function(i)
    with(tasks[i, ], tune_task(dataset, hazard, learner)),
    mc.cores = cores, mc.preschedule = FALSE)
  tuned <- list()
  for (o in out) tuned[[o$dataset]][[o$hazard]][[o$learner]] <- o$params
  dir.create("results", showWarnings = FALSE)
  saveRDS(tuned, "results/tuned_armA.rds")
  log <- do.call(rbind, lapply(out, function(o) data.frame(
    dataset = o$dataset, hazard = o$hazard, learner = o$learner,
    params = if (is.null(o$params)) NA else paste(names(o$params), unlist(o$params), sep = "=", collapse = "; "),
    ibs = if (is.null(o$ibs)) NA else o$ibs, secs = round(o$secs), error = o$error)))
  write.csv(log, "results/tuned_armA.csv", row.names = FALSE)
  print(log, row.names = FALSE)
}

# install-deps.R -- install the exact package versions used for the run.
# CRAN packages first, then the ielbadisy GitHub packages (main branch) in
# dependency order. Run once:  Rscript install-deps.R

repos <- "https://cloud.r-project.org"
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = repos)

cran <- c("survival", "CoxBoost", "flexsurv", "ranger", "randomForestSRC", "aorsf",
          "partykit", "party", "mboost", "survivalsvm", "glmnet", "pec", "aftgee",
          "nnls", "pracma", "gower", "BART", "timereg", "xgboost", "rsample",
          "rpart", "ggplot2", "glue", "cli", "rlang")
miss <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss, repos = repos)

# version each GitHub package had when the run started (2026-11-15)
gh <- c(basetable = "1.4.8", functionals = "0.6.0", unsurv = "0.7.3", tvrmst = "0.0.9",
        densemlp = "0.9.3", fastgbm = "0.7.3", rpsurv = "0.8.4", nbsurv = "0.5.4",
        mcstatsim = "0.9.6", survalis = "2.4.2")
for (p in names(gh)) {
  have <- tryCatch(as.character(packageVersion(p)), error = function(e) "")
  if (have != gh[[p]]) remotes::install_github(paste0("ielbadisy/", p), upgrade = "never",
                                               dependencies = TRUE, repos = repos)
}

# check
for (p in names(gh)) {
  v <- as.character(packageVersion(p))
  cat(sprintf("%-12s %-8s %s\n", p, v, if (v == gh[[p]]) "ok" else paste("EXPECTED", gh[[p]])))
}

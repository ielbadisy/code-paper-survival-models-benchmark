# install-deps.R -- install the packages used for the published run.
# CRAN packages first, then the author's GitHub packages at the versions used
# (2026-11-17). Run once from the bundle root:  Rscript install-deps.R

repos <- "https://cloud.r-project.org"
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = repos)
cran <- c("survival", "CoxBoost", "glmnet", "randomForestSRC", "aorsf", "mboost", "BART", "xgboost",
          "ggplot2", "ggrepel", "patchwork")
miss <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss, repos = repos)

GH <- c(basetable = "1.4.9", functionals = "0.6.0", unsurv = "0.7.3", tvrmst = "0.0.9",
        densemlp = "0.9.3", fastgbm = "0.7.3", rpsurv = "0.8.4", nbsurv = "0.6.0",
        mcstatsim = "0.9.6", survalis = "2.4.3")
for (p in names(GH)) if (!requireNamespace(p, quietly = TRUE) || packageVersion(p) != GH[[p]])
  remotes::install_github(paste0("ielbadisy/", p), upgrade = "never", dependencies = TRUE)
for (p in names(GH)) cat(sprintf("%-12s %-8s %s\n", p, packageVersion(p),
                                 if (packageVersion(p) == GH[[p]]) "ok" else paste("expected", GH[[p]])))

# Code and results: survival learners in randomised trial data

Code and per-replicate results for "Beyond Cox versus random forests: a neutral comparison of
survival learners for patient-specific prediction and treatment benefit in randomised trial data"
(Imad El Badisy). Every table and figure of the paper and of Supporting Information S1 is produced
from the files below.
The same files are submitted as Supporting Information S2 of the paper.

## Layout

- `R/`: data-generating mechanisms (`dgp.R`), performance and benefit measures (`metrics.R`) and
  shared helpers for the tables and figures (`res-utils.R`).
- `data-ref/adv_prostate_ca/`: the VACURG advanced prostate cancer trial (Byar and Green 1980).
  The PBC trial is taken from the R package survival.
- `armA-rct/`: the simulation study.
  - `dgp-check.R`: large-sample checks of the data-generating mechanisms (censoring share, true
    survival against latent times, Kaplan-Meier bias). Output `results/dgp_check.csv`.
  - `common.R`: design, learners and scoring shared by tuning and evaluation.
  - `tune-run.R`: tuning on replicates 1 to 10 (seed 1027). Output `results/tuned_v2.rds` and
    `results/tune_raw_v2.rds`.
  - `sim-run.R`: evaluation on 90 replicates per scenario (seed 2028). Output
    `results/sim_raw_v2_90reps.rds` (one row per scenario, replicate, learner and measure).
  - `sim-res.R`: means, Monte Carlo standard errors, ranks, failure rates and timing.
  - `explore-res.R`: paired differences from Cox PH, metric selection, prediction against benefit,
    policy value and true benefit spread (`results/analysis_v2_90reps/`).
  - `fig-res.R`: Tables 2 and 3 and Figures 1 to 5 (`results/figures_v2_90reps/`) on the 24 core
    scenarios (n = 200, 30% censoring), plus full-design tables for S1.
  - `fig-extra.R`: heatmap and nested loop plot of Supporting Information S1.
  - `decomp-run.R` and `decomp-res.R`: decomposition of benefit error into prediction error and
    cross-arm correlation (30 replicates of the 24 base scenarios with n = 200, seed 3029).
  - `dca-run.R` and `dca-res.R`: treatment decision curves against the truth (30 replicates of the
    same 24 scenarios, seed 4031) and their figure in Supporting Information S1.
- `armB-real/`: the ACTG 175 example. `actg175-run.R` tunes and evaluates on 2:1 splits,
  `actg175-res.R` makes Table 4 and Figures 6 to 8, and `actg175-null.R` scores a Kaplan-Meier curve
  without covariates on the same splits (last row of Table 4). `actg175-dca.R` computes the treatment and
  prognostic decision curves on the same splits, and `actg175-dca-res.R` draws Figure 9 and its
  supplementary figure. The data are taken from the R package BART.
- `paper-numbers.R`: every number quoted in the text of the paper and of Supporting Information S1,
  recomputed from the stored results (`results/paper_numbers.csv`). Run from the repository root.

## Reproduce

Install the packages once:

```bash
Rscript install-deps.R
```

Remake all summaries, tables and figures from the stored raw results (a few minutes):

```bash
cd armA-rct
Rscript sim-res.R v2_90reps
Rscript explore-res.R v2_90reps
Rscript fig-res.R v2_90reps
Rscript fig-extra.R v2_90reps
cd ../armB-real
Rscript actg175-res.R 50
Rscript actg175-null.R
Rscript actg175-dca-res.R
cd ../armA-rct
Rscript dca-res.R
cd ..
Rscript paper-numbers.R
```

Rerun the experiments from scratch. Pin the BLAS threads to 1 so that parallel workers do not
oversubscribe the cores:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd armA-rct
Rscript tune-run.R <cores>          # tuning, under 10 CPU hours
Rscript sim-run.R 90 <cores>        # evaluation, about 60 CPU hours
cd ../armB-real
Rscript actg175-run.R tune <cores>
Rscript actg175-run.R eval <cores> 50
Rscript actg175-dca.R <cores>
cd ../armA-rct
Rscript decomp-run.R 30 <cores>
Rscript dca-run.R 30 <cores>
Rscript dgp-check.R
```

All runs use fixed seeds through mcstatsim and give the same results on the same package versions.
The analyses used R 4.5.1 on Linux.

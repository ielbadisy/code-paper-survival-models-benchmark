# code-paper-survival-models-benchmark

Simulation code for a trial-shaped benchmark of survival learners (survalis),
with prediction metrics, calibration against the true survival function, and
individual treatment-benefit metrics. The DGPs are built on two reference
trials: PBC (`survival::pbc`) and VACURG advanced prostate cancer
(`data-ref/adv_prostate_ca/`).

## Layout

- `R/` DGPs (`dgp-rct.R`, `dgp-rct-stress.R`), metrics (`sim-core.R`, `benefit.R`)
- `armA-rct/tune-run.R` tunes each learner once per dataset x hazard (output already in `results/tuned_armA.rds`)
- `armA-rct/sim-run.R` full run: 156 cells x 200 replicates x 15 learners
- `armA-rct/sim-res.R` summaries (means, MCSE, ranks, failure rates, timing)
- `armA-rct/example-res.R` real-data illustration

## Launch (copy-paste)

On the workstation:

```bash
git clone git@github.com:ielbadisy/code-paper-survival-models-benchmark.git
cd code-paper-survival-models-benchmark
Rscript install-deps.R     # once; each line should say "ok"
./run.sh 30                # number of cores; default is all minus 2
tail -c 300 armA-rct/results/sim-run.log
```

`run.sh` pins `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS`
to 1 before R starts. Do not launch `sim-run.R` without these.

Cost: about 24 CPU seconds per cell x replicate, so about 210 CPU hours in
total for 31 200 cell-replicates. The run checkpoints every 5 replicates to
`armA-rct/results/checkpoint_200reps.rds`. If it stops, run `./run.sh` again
and it resumes. The committed checkpoint already holds replicates 1 to 70
(about 35%), so a fresh clone starts at replicate 71.

## When it finishes

```bash
cd armA-rct
Rscript sim-res.R 200reps
Rscript example-res.R
```

Outputs land in `armA-rct/results/` (`sim_raw_200reps.{rds,csv}`,
`sim_summary_200reps.csv`, `sim_mean_rank_200reps.csv`,
`sim_fail_rates_200reps.csv`, `timing_200reps.csv`). Report
`benefit_cal_slope` by its median (heavy-tailed).

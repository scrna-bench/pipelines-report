# pipelines-report

REPORT stage of the scRNA pipelines benchmark: an omnibenchmark (api 0.8)
global gather over every METRIC node. `report.R` concatenates the per-dataset
`metrics_tsv` and `timings_tsv` tables into `metrics.tsv` and `timings.tsv`
and renders `report.Rmd` to a self-contained `plots.html`.

`envs/r/` is an [rv](https://github.com/A2-ai/rv) project; CI builds it into
`ghcr.io/scrna-bench/pipelines-report:r-<sha>`.

#!/usr/bin/env Rscript
# REPORT gather (omnibenchmark api 0.8): concatenates every METRIC node's
# metrics_tsv and timings_tsv (one node per dataset) and renders report.Rmd.
#
#   --from-files <lineage.json>   the METRIC nodes (gather manifest)
#   --output report_metrics_tsv=<path> --output report_timings_tsv=<path>
#   --output report_html=<path>

here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "src", "cli.R"))

args <- parse_cli()
require_args(args, "from-files")
out <- args$outputs
for (id in c("report_metrics_tsv", "report_timings_tsv", "report_html")) {
  if (is.null(out[[id]])) stop("missing --output ", id, "=PATH", call. = FALSE)
}

members <- jsonlite::fromJSON(args[["from-files"]], simplifyVector = FALSE)$members
read_all <- function(id) {
  tables <- lapply(members, function(m) {
    utils::read.delim(m$files[[id]], check.names = FALSE, stringsAsFactors = FALSE)
  })
  dplyr::bind_rows(tables)
}
metrics <- read_all("metrics_tsv")
timings <- read_all("timings_tsv")

write.table(metrics, out$report_metrics_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(timings, out$report_timings_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

rmarkdown::render(
  file.path(here, "report.Rmd"),
  output_file = basename(out$report_html),
  output_dir = dirname(out$report_html),
  intermediates_dir = tempdir(),
  params = list(metrics = metrics, timings = timings),
  quiet = TRUE,
  envir = new.env(parent = globalenv())
)

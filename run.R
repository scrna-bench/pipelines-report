#!/usr/bin/env Rscript

library(argparse)
library(jsonlite)

parser <- ArgumentParser(description = "Benchmarking entrypoint")

parser$add_argument(
  "--output_dir", "-o",
  dest = "output_dir", type = "character",
  help = "output directory where files will be saved",
  default = getwd(), required = TRUE
)
parser$add_argument(
  "--name", "-n",
  dest = "name", type = "character",
  help = "name of the module",
  required = TRUE
)
parser$add_argument(
  "--metrics.json",
  dest = "metrics_paths", type = "character", nargs = "+",
  help = "cluster tsv path", required = TRUE
)

args <- parser$parse_args()

cargs <- commandArgs(trailingOnly = FALSE)
m <- grep("--file=", cargs)
run_dir <- dirname(gsub("--file=", "", cargs[[m]]))

metrics_path <- file.path(args$output_dir, "metrics.tsv")
timings_path <- file.path(args$output_dir, "timings.tsv")
report_path <- file.path(args$output_dir, "plots.html")

split_path <- function(p) strsplit(p, "/", fixed = TRUE)[[1]]

# function to extract dataset name, method name
# and clustering resolution from path
extract_run_info <- function(p) {
  parts <- split_path(p)

  i_ds <- match("datasets", parts)
  dataset_dir <- paste(parts[1:(i_ds + 2)], collapse = "/")
  cfg <- fromJSON(file.path(dataset_dir, "parameters.json"))
  dataset_name <- cfg$dataset_name

  i_methods <- match("methods", parts)
  method_dir <- paste(parts[1:(i_methods + 2)], collapse = "/")
  cfg <- fromJSON(file.path(method_dir, "parameters.json"))
  method_name <- cfg$method_name
  resolution <- cfg$resolution
  filtering <- cfg$filter
  n_comp <- cfg$n_comp
  n_neig <- cfg$n_neig

  data.frame(
    dataset = dataset_name, method = method_name,
    resolution = resolution, filtering = filtering,
    n_comp = n_comp, n_neig = n_neig
  )
}

metrics <- vector("list", length(args$metrics_paths))
timings <- vector("list", length(args$metrics_paths))

# build up ARI and timing rows in lists
for (i in seq_along(args$metrics_paths)) {
  p <- args$metrics_paths[i]
  x <- fromJSON(p)

  meta <- extract_run_info(p)

  row <- list()
  row$n_clusters_leiden <- x$n_clusters$leiden
  row$n_clusters_louvain <- x$n_clusters$louvain
  row$dropped_cells <- x$dropped_cells

  for (metric_name in names(x$agreement)) {
    for (algo in names(x$agreement[[metric_name]])) {
      key <- paste("agree", metric_name, algo, sep = "_")
      row[[key]] <- x$agreement[[metric_name]][[algo]]
    }
  }
  for (metric_name in names(x$structure)) {
    for (algo in names(x$structure[[metric_name]])) {
      key <- paste("struc", metric_name, algo, sep = "_")
      row[[key]] <- x$structure[[metric_name]][[algo]]
    }
  }

  metrics[[i]] <- cbind(meta, as.data.frame(row, check.names = FALSE))

  tdf <- as.data.frame(x$timings)
  timings[[i]] <- cbind(meta, tdf)
}

# rbind lists into dataframes
metrics <- do.call(rbind, metrics)
timings <- do.call(rbind, timings)

write.table(
  metrics, metrics_path,
  sep = "\t", quote = F, row.names = F
)
write.table(
  timings, timings_path,
  sep = "\t", quote = F, row.names = F
)

rmarkdown::render(
  file.path(run_dir, "plots.Rmd"),
  output_file = "plots.html",
  output_dir = args$output_dir,
  params = list(
    metrics = metrics,
    timings = timings
  ),
  quiet = TRUE,
  envir = new.env(parent = globalenv())
)

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

ari_path <- file.path(args$output_dir, "ari.tsv")
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

  data.frame(
    dataset = dataset_name, method = method_name,
    resolution = resolution, filtering = filtering
  )
}

ari <- vector("list", length(args$metrics_paths))
timings <- vector("list", length(args$metrics_paths))

# build up ARI and timing rows in lists
for (i in seq_along(args$metrics_paths)) {
  p <- args$metrics_paths[i]
  x <- fromJSON(p)

  meta <- extract_run_info(p)

  adf <- as.data.frame(x$ari)
  adf$n_clusters_leiden <- x$n_clusters$leiden
  adf$n_clusters_louvain <- x$n_clusters$louvain
  ari[[i]] <- cbind(meta, adf)

  tdf <- as.data.frame(x$timings)
  timings[[i]] <- cbind(meta, tdf)
}

# rbind lists into dataframes
ari <- do.call(rbind, ari)
timings <- do.call(rbind, timings)

write.table(
  ari, ari_path,
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
    ari = ari,
    timings = timings
  ),
  quiet = TRUE,
  envir = new.env(parent = globalenv())
)

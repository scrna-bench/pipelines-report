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

ari_path <- file.path(args$output_dir, "ari.tsv")
timings_path <- file.path(args$output_dir, "timings.tsv")

split_path <- function(p) strsplit(p, "/", fixed = TRUE)[[1]]

# function to extract dataset and method name from path
extract_ds_method <- function(p) {
  parts <- split_path(p)

  i_ds <- match("datasets", parts)
  dataset_dir <- paste(parts[1:(i_ds + 2)], collapse = "/")
  cfg <- fromJSON(file.path(dataset_dir, "parameters.json"))
  dataset_name <- cfg$dataset_name

  i_methods <- match("methods", parts)
  method_dir <- paste(parts[1:(i_methods + 2)], collapse = "/")
  cfg <- fromJSON(file.path(method_dir, "parameters.json"))
  method_name <- cfg$method_name

  list(dataset = dataset_name, method = method_name)
}

ari <- vector("list", length(args$metrics_paths))
timings <- vector("list", length(args$metrics_paths))

# build up ARI and timing rows in lists
for (i in seq_along(args$metrics_paths)) {
  p <- args$metrics_paths[i]
  x <- fromJSON(p)

  meta <- extract_ds_method(p)

  adf <- as.data.frame(x$ari)
  ari[[i]] <- cbind(
    data.frame(dataset = meta$dataset, method = meta$method),
    adf
  )

  tdf <- as.data.frame(x$timings)
  timings[[i]] <- cbind(
    data.frame(dataset = meta$dataset, method = meta$method),
    tdf
  )
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

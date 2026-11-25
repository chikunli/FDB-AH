entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(entry)), "lib", "init.R"))
sim <- new.env(parent = baseenv())
sys.source(file.path(repro_dir, "lib", "simulation.R"), sim)

design <- subset(sim$simulation_design(), study == "clean")
raw <- sim$simulate(design, repro_map)
methods <- c("FDB-AH", "FDB-FS", "FDB-HF", "FDB-75")
out <- expand.grid(method = methods, setting = design$setting, stringsAsFactors = FALSE)
out <- out[c("setting", "method")]
for (metric in sim$metric_names) {
  out[[metric]] <- vapply(seq_len(nrow(out)), function(i) {
    id <- design$id[match(out$setting[i], design$setting)]
    cell <- sim$table_entries(raw, id, out$method[i], metric)
    sprintf("%s (%s)", cell["mean"], cell["sd"])
  }, character(1L))
}
utils::write.csv(out, file.path(output_dir, "Table1.csv"), row.names = FALSE)

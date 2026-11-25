entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(entry)), "lib", "init.R"))
sim <- new.env(parent = baseenv())
sys.source(file.path(repro_dir, "lib", "simulation.R"), sim)

design <- subset(sim$simulation_design(), study == "fixed" & eps == 0.10)
raw <- sim$simulate(design, repro_map)
types <- names(sim$type_labels)
out <- expand.grid(statistic = c("mean", "sd"), metric = sim$metric_names,
                   type = types, stringsAsFactors = FALSE)
out <- out[c("type", "metric", "statistic")]
for (setting in c("A", "B", "C")) {
  for (method in c("FDB-AH", "FDB-FS", "FDB-75")) {
    column <- paste(setting, method, sep = "_")
    out[[column]] <- vapply(seq_len(nrow(out)), function(i) {
      id <- design$id[design$setting == setting & design$type == out$type[i]]
      sim$table_entries(raw, id, method, out$metric[i])[out$statistic[i]]
    }, character(1L), USE.NAMES = FALSE)
  }
}
out$type <- unname(sim$type_labels[out$type])
utils::write.csv(out, file.path(output_dir, "Table2.csv"), row.names = FALSE)

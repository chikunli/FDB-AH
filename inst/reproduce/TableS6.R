entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(entry)), "lib", "init.R"))
source(file.path(repro_dir, "lib", "real.R"))

tables <- lapply(c("normal", "rotor"), function(record) {
  runs <- run_record(record)
  data.frame(record = record, record_means(runs$scores, c("h", "FPR", "FNR", "seconds")))
})
utils::write.csv(do.call(rbind, tables), file.path(output_dir, "TableS6.csv"), row.names = FALSE)

entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(entry)), "lib", "init.R"))
source(file.path(repro_dir, "lib", "real.R"))

skin <- run_record("skin")
columns <- c("h", "initial_FPR", "initial_FNR", "FPR", "FNR", "seconds")
utils::write.csv(record_means(skin$scores, columns),
                 file.path(output_dir, "Table4.csv"), row.names = FALSE)

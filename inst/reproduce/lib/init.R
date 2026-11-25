repro_dir <- dirname(normalizePath(entry, mustWork = TRUE))
output_dir <- Sys.getenv("FDBAH_OUTPUT", unset = file.path(getwd(), "output"))
if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE))
    stop("Cannot create the output directory: ", output_dir)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
if (!requireNamespace("FDBAH", quietly = TRUE)) stop("Install the FDBAH package first.")

repro_map <- function(X, FUN) {
    workers <- as.numeric(Sys.getenv("FDBAH_WORKERS", "1"))
    if (length(workers) != 1L || !is.finite(workers) ||
        workers < 1 || workers != floor(workers))
        stop("FDBAH_WORKERS must be a positive integer.")
    workers <- min(workers, length(X))
    if (workers <= 1L) return(lapply(X, FUN))
    cl <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cl))
    helper <- if (basename(entry) %in% c("Table1.R", "Table2.R", "Figure3.R"))
        "simulation.R" else "real.R"
    parallel::clusterCall(cl, function(libs, scripts, out, helper) {
        .libPaths(libs)
        assign("repro_dir", scripts, envir = .GlobalEnv)
        assign("output_dir", out, envir = .GlobalEnv)
        sys.source(file.path(scripts, "lib", helper), envir = .GlobalEnv)
        NULL
    }, .libPaths(), repro_dir, output_dir, helper)
    parallel::parLapplyLB(cl, X, FUN)
}

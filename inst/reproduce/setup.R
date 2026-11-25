entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
scripts <- dirname(normalizePath(entry, mustWork = TRUE))
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% "--comparators")) stop("Usage: Rscript setup.R [--comparators]")
comparators <- "--comparators" %in% args
cran <- c("ddalpha", "robustbase", "MASS", "digest", "data.table", "ggplot2")
if (comparators) cran <- c(cran, "Rcpp", "RcppArmadillo", "BH")
missing <- cran[!vapply(cran, requireNamespace, logical(1L), quietly = TRUE)]
if (requireNamespace("ggplot2", quietly = TRUE) && utils::packageVersion("ggplot2") < "3.5.0")
    missing <- union(missing, "ggplot2")
if (length(missing)) utils::install.packages(missing, repos = "https://cloud.r-project.org")
if (!all(vapply(cran, requireNamespace, logical(1L), quietly = TRUE)))
    stop("Some required R packages are unavailable.")

install_source <- function(path) {
    withCallingHandlers(
        utils::install.packages(path, repos = NULL, type = "source", INSTALL_opts = "--preclean"),
        warning = function(w) {
            if (grepl("non-zero exit status", conditionMessage(w), fixed = TRUE))
                stop("Installation failed for ", basename(path), call. = FALSE)
        }
    )
}

install_fsdac <- function() {
    commit <- "85ede5feeb450ac0329ae0c1e420d9430aee8980"
    zip <- Sys.getenv("FDBAH_FSDA_ZIP")
    if (!nzchar(zip)) {
        zip <- tempfile(fileext = ".zip")
        utils::download.file(paste0("https://codeload.github.com/UniprJRC/",
                                    "FSDA-MATLAB_Coder/zip/", commit), zip, mode = "wb")
    }
    members <- utils::unzip(zip, list = TRUE)$Name
    prefix <- paste0("FSDA-MATLAB_Coder-", commit, "/R-package/")
    part <- members[startsWith(members, prefix)]
    if (!length(part) || any(grepl("(^|/)\\.\\.(/|$)", part)))
        stop("Expected the FSDA archive for commit ", commit)
    work <- tempfile("fsdac-")
    dir.create(work)
    utils::unzip(zip, files = part, exdir = work)
    pkg <- file.path(work, prefix)
    unused <- c("src/_fsminvmmd.c", "src/_fsrinvmdr.c", "R/FSMinvmmd.R", "R/FSRinvmdr.R")
    if (!all(file.remove(file.path(pkg, unused)))) stop("Could not prepare the FSDA build.")
    reg <- file.path(pkg, "src", "_fsdac_init.c")
    lines <- readLines(reg, warn = FALSE)
    omit <- grepl("r_fsminvmmd|r_fsrinvmdr", lines)
    stopifnot(sum(omit) == 4L)
    writeLines(lines[!omit], reg)
    ns <- file.path(pkg, "NAMESPACE")
    exports <- readLines(ns, warn = FALSE)
    writeLines(gsub("FSMinvmmd, |FSRinvmdr, ", "", exports), ns)
    install_source(pkg)
}

if (comparators) {
    stable <- file.path(scripts, "vendor", "StableMCD_1.0.tar.gz")
    hash <- digest::digest(file = stable, algo = "sha256", serialize = FALSE)
    stopifnot(hash == "1001e153ff57e075cb197b88b02be44b2ec6a9ed09960112d64378d1ed48aebb")
    install_source(stable)
    install_fsdac()
    stopifnot(as.character(utils::packageVersion("StableMCD")) == "1.0",
              utils::packageVersion("fsdac") == "0.6-1")
}

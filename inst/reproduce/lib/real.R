read_record <- function(record) {
  work <- tempfile("fdbah-raw-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  python <- Sys.getenv("FDBAH_PYTHON")
  if (!nzchar(python)) {
    found <- Sys.which(c("python3", "python"))
    python <- unname(found[nzchar(found)])[1L]
  }
  if (!length(python) || is.na(python) || !nzchar(Sys.which(python))) {
    stop("Set FDBAH_PYTHON to a Python 3 interpreter with numpy and pandas.")
  }
  csv <- file.path(work, "sample.csv")
  args <- c(file.path(repro_dir, "lib", "raw.py"), record, csv)
  status <- system2(python, shQuote(args))
  if (status != 0L) {
    stop("Raw-data preparation failed. Check Python (numpy and pandas) and the source files.")
  }
  data <- data.table::fread(csv, data.table = FALSE, nThread = 1L)
  label <- data$label
  x <- as.matrix(data[setdiff(names(data), "label")])
  storage.mode(x) <- "double"
  expected <- switch(record, normal = c(1500, 8, 0), rotor = c(1090, 8, 347),
                     skin = c(2500, 3, 491))
  stopifnot(all(is.finite(x)), all(label %in% 0:1),
            all(c(dim(x), sum(label)) == expected))
  center <- apply(x, 2L, stats::median)
  scale <- apply(x, 2L, robustbase::Qn)
  fallback <- !is.finite(scale) | scale <= 0
  if (any(fallback)) scale[fallback] <- apply(x[, fallback, drop = FALSE], 2L, stats::sd)
  stopifnot(all(is.finite(scale) & scale > 0))
  list(x = sweep(sweep(x, 2L, center, "-"), 2L, scale, "/"), label = label)
}

real_methods <- c("FDB-AH", "FDB-FS", "FDB-HF", "FDB-75", "FSDA", "Stable MCD")

native_fit <- function(x, method, seed) {
  set.seed(seed)
  start <- Sys.time()
  if (method == "FSDA") {
    fit <- fsdac::FSM(x, msg = FALSE)
    out <- as.integer(fit$outliers)
    stopifnot(all(out %in% seq_len(nrow(x))))
    initial <- NULL
    retained <- setdiff(seq_len(nrow(x)), out)
  } else {
    alpha <- seq(0.5, 0.975, 0.025)
    selection <- StableMCD::bootstrap_mcd(x, alpha, B = 50L, classifier = "MD", sd_ratio = 3)
    stopifnot(length(selection$best_alpha) > 0L, all(is.finite(selection$iim)))
    fit <- StableMCD::mcd(x, min(selection$best_alpha), verbose = FALSE, reweighting = FALSE)
    initial <- retained <- fit$index
  }
  list(initial = initial, retained = retained, elapsed = as.numeric(difftime(Sys.time(), start, units = "secs")))
}

record_scores <- function(fit, method, label, seed) {
  error_rates <- function(kept) {
    if (is.null(kept)) return(c(NA_real_, NA_real_))
    stopifnot(!anyDuplicated(kept), all(kept %in% seq_along(label)))
    inside <- seq_along(label) %in% kept
    c(100 * mean(!inside[label == 0L]),
      if (any(label == 1L)) 100 * mean(inside[label == 1L]) else NA_real_)
  }
  initial <- error_rates(fit$initial)
  final <- error_rates(fit$retained)
  data.frame(method = method, seed = seed,
             h = if (is.null(fit$initial)) NA_integer_ else length(fit$initial),
             initial_FPR = initial[1], initial_FNR = initial[2],
             FPR = final[1], FNR = final[2],
             seconds = if (is.null(fit$elapsed)) NA_real_ else fit$elapsed)
}

run_record <- function(record, figure = FALSE) {
  data <- read_record(record)
  seeds <- if (record == "skin") 101:200 else 1:100
  runs <- repro_map(seeds, function(seed) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(seed)
    if (figure) {
      methods <- c("FDB-AH", "FDB-HF", "FDB-75")
      fits <- lapply(methods, function(m) FDBAH::fdb_fit(data$x, method = m, seed = seed))
      names(fits) <- methods
    } else {
      if (!requireNamespace("StableMCD", quietly = TRUE) ||
          !requireNamespace("fsdac", quietly = TRUE)) stop("Install StableMCD and fsdac first.")
      fits <- FDBAH::fdb_four(data$x, seed = seed, timed = TRUE)
      for (m in c("FSDA", "Stable MCD")) fits[[m]] <- native_fit(data$x, m, seed)
    }
    scores <- lapply(names(fits), function(m) record_scores(fits[[m]], m, data$label, seed))
    list(scores = do.call(rbind, scores),
         path = if (figure && seed == 1L) fits[["FDB-AH"]]$path else NULL)
  })
  rows <- do.call(rbind, lapply(runs, "[[", "scores"))
  stopifnot(all(table(rows$method) == length(seeds)))
  list(scores = rows, path = runs[[1L]]$path)
}

record_means <- function(rows, metrics) {
  ans <- do.call(rbind, lapply(real_methods, function(m) {
    trials <- rows[rows$method == m, metrics, drop = FALSE]
    stopifnot(nrow(trials) == 100L)
    data.frame(method = m, t(vapply(trials, mean, 0.0)), check.names = FALSE)
  }))
  for (metric in metrics) {
    x <- ans[[metric]]
    ans[[metric]] <- ifelse(is.na(x), "---", sprintf("%.3f", x))
  }
  ans
}

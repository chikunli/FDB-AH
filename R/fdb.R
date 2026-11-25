fdb_methods <- c("FDB-AH", "FDB-FS", "FDB-HF", "FDB-75")

prepare_sample <- function(x, seed, directions, window) {
    x <- as.matrix(x)
    if (!is.numeric(x) || ncol(x) < 1L)
        stop("x must be a numeric matrix with at least one column.")
    rows <- which(stats::complete.cases(x))
    if (length(rows) != nrow(x)) {
        warning("Omitting rows with missing coordinates.")
        x <- x[rows, , drop = FALSE]
    }
    n <- nrow(x)
    p <- ncol(x)
    if (n < p + 2L || any(!is.finite(x)))
        stop("At least p + 2 complete rows and finite coordinates are required.")
    if (!is.numeric(window) || length(window) != 1L || !is.finite(window) ||
        window <= 0 || window >= 1)
        stop("window must be between zero and one.")
    if (is.null(directions)) directions <- max(1000L, 10L * p)
    positive_integer <- function(value) {
        is.numeric(value) && length(value) == 1L && is.finite(value) &&
            value >= 1 && value <= .Machine$integer.max && value == floor(value)
    }
    if (!positive_integer(seed) || !positive_integer(directions))
        stop("seed and directions must be positive integers.")
    depth <- as.numeric(ddalpha::depth.projection(
        x, x, num.directions = as.integer(directions), seed = as.integer(seed)))
    if (length(depth) != n || any(!is.finite(depth)))
        stop("The projection depth is undefined for this sample.")
    list(x = x, order = order(depth, decreasing = TRUE), rows = rows,
         n = n, p = p, half = floor((n + p + 1) / 2),
         radius = max(1L, floor(window * n)))
}

subset_model <- function(sample, h, fixed_simulation = FALSE) {
    x <- sample$x
    first <- sample$order[seq_len(h)]
    subset <- x[first, , drop = FALSE]
    center <- colMeans(subset)
    consistency <- robustbase::.MCDcons(sample$p, h / sample$n)
    finite_sample <- robustbase::.MCDcnp2(sample$p, sample$n, h / sample$n)
    scatter <- stats::cov(subset)
    if (fixed_simulation) {
        scatter <- (consistency * finite_sample) * scatter
    } else {
        scatter <- scatter * consistency * finite_sample
    }
    distances <- tryCatch(stats::mahalanobis(x, center, scatter),
                          error = function(e) rep(NA_real_, sample$n))
    list(center = center, cov = scatter, distance = distances)
}

fit_subset <- function(sample, h, path = NULL, fixed_simulation = FALSE) {
    model <- subset_model(sample, h, fixed_simulation)
    distance <- model$distance
    cutoff <- stats::qchisq(0.975, sample$p)
    if (fixed_simulation) {
        keep <- which(is.finite(distance) & distance < cutoff)
    } else {
        keep <- which(sqrt(distance) <= sqrt(cutoff))
    }
    center <- model$center
    scatter <- model$cov
    if (length(keep) <= sample$p) {
        warning("Too few observations retained; returning the initial estimates.")
    } else if (fixed_simulation) {
        final <- sample$x[keep, , drop = FALSE]
        center <- colMeans(final)
        scatter <- stats::cov(final)
    } else {
        weights <- numeric(sample$n)
        weights[keep] <- 1
        final <- stats::cov.wt(sample$x, wt = weights, center = TRUE)
        center <- final$center
        scatter <- final$cov
    }
    list(h = h, initial = sample$rows[sample$order[seq_len(h)]],
         retained = sample$rows[keep], center = center, cov = scatter,
         raw_center = model$center, raw_cov = model$cov,
         raw_distance = distance, path = path)
}

fit_method <- function(sample, method) {
    switch(method,
        "FDB-AH" = adaptive_fits(sample, both = FALSE, cache = FALSE)[[1L]],
        "FDB-FS" = {
            path <- direct_path(sample, stop_at_crossing = TRUE)
            fit_subset(sample, path$h[signal_index(path)], path)
        },
        "FDB-HF" = fit_subset(sample, sample$half),
        "FDB-75" = fit_subset(sample, max(sample$p + 1L, floor(0.75 * sample$n)))
    )
}

fdb_fit <- function(x, method = "FDB-AH", seed = 1L, directions = NULL,
                    window = 0.08) {
    method <- match.arg(method, fdb_methods)
    sample <- prepare_sample(x, seed, directions, window)
    fit_method(sample, method)
}

fdb_four <- function(x, seed = 1L, timed = FALSE, directions = NULL,
                     window = 0.08) {
    if (!is.logical(timed) || length(timed) != 1L || is.na(timed))
        stop("timed must be TRUE or FALSE.")
    start <- Sys.time()
    sample <- prepare_sample(x, seed, directions, window)
    depth_seconds <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    if (timed) {
        fits <- lapply(fdb_methods, function(method) {
            start <- Sys.time()
            fit <- fit_method(sample, method)
            fit$elapsed <- depth_seconds +
                as.numeric(difftime(Sys.time(), start, units = "secs"))
            fit
        })
    } else {
        adaptive <- adaptive_fits(sample, both = TRUE, cache = TRUE)
        fixed_sizes <- c(sample$half, max(sample$p + 1L, floor(0.75 * sample$n)))
        fixed <- lapply(fixed_sizes, function(h) {
            fit_subset(sample, h, fixed_simulation = TRUE)
        })
        fits <- c(adaptive, fixed)
    }
    names(fits) <- fdb_methods
    fits
}

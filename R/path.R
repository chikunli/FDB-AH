signal_limit <- function(h, n, p) {
    q <- stats::qf(1 - 0.999, 2 * (n - h), 2 * (h + 1))
    u <- (h + 1) / (h + 1 + (n - h) * q)
    factor <- p * ((h + 1) / h) * ((h - 1) / (h - p))
    sqrt(factor * stats::qf(u, p, h - p))^2
}

constants_cache <- new.env(parent = emptyenv())

path_constants <- function(sample, cache) {
    n <- sample$n
    p <- sample$p
    if (cache && identical(c(n, p), constants_cache$dimensions))
        return(constants_cache$value)
    h <- seq.int(sample$half, n - 1L)
    c1 <- vapply(h, function(size) robustbase::.MCDcons(p, size / n), numeric(1))
    c2 <- vapply(h, function(size) robustbase::.MCDcnp2(p, n, size / n), numeric(1))
    constants <- data.frame(h = h, c1 = c1, c2 = c2,
                            scale = c1 * c2, threshold_sq = signal_limit(h, n, p))
    if (cache) {
        constants_cache$dimensions <- c(n, p)
        constants_cache$value <- constants
    }
    constants
}

rank_distance <- function(distance, h) {
    finite <- distance[is.finite(distance)]
    if (length(finite) <= h) return(NA_real_)
    sort(finite, partial = h + 1L)[h + 1L]
}

direct_path <- function(sample, stop_at_crossing = FALSE) {
    sizes <- seq.int(sample$half, sample$n - 1L)
    limits <- signal_limit(sizes, sample$n, sample$p)
    distance <- rep(NA_real_, length(sizes))
    for (j in seq_along(sizes)) {
        h <- sizes[j]
        distance[j] <- rank_distance(subset_model(sample, h)$distance, h)
        if (stop_at_crossing && is.finite(distance[j]) && distance[j] > limits[j])
            break
    }
    used <- seq_len(j)
    data.frame(h = sizes[used], d_next_sq = distance[used],
               threshold_sq = limits[used])
}

signal_index <- function(path) {
    crossed <- which(is.finite(path$d_next_sq) & path$d_next_sq > path$threshold_sq)
    if (length(crossed)) return(crossed[1L])
    last_valid(path$d_next_sq)
}

last_valid <- function(distance) {
    valid <- which(is.finite(distance))
    if (!length(valid)) stop("All candidate subsets gave undefined ordered distances.")
    valid[length(valid)]
}

window_choice <- function(distance, radius) {
    best <- NA_integer_
    near_tie <- FALSE
    for (j in seq_along(distance)) {
        window <- seq.int(max(1L, j - radius), min(length(distance), j + radius))
        neighbours <- distance[window[window != j]]
        if (!is.finite(distance[j]) || any(!is.finite(neighbours))) next
        margin <- if (length(neighbours)) distance[j] - max(neighbours) else Inf
        if (is.na(best) && margin >= 0) best <- j
        if (abs(margin) <= 1e-6 * max(1, abs(distance[j]))) near_tie <- TRUE
    }
    if (is.na(best)) best <- last_valid(distance)
    list(index = best, near_tie = near_tie)
}

recursive_path <- function(sample, cache) {
    x <- sample$x
    constants <- path_constants(sample, cache)
    sizes <- constants$h
    distance <- numeric(length(sizes))
    refresh_error <- 0
    for (j in seq_along(sizes)) {
        h <- sizes[j]
        scale <- constants$scale[j]
        if ((j - 1L) %% 64L == 0L) {
            initial <- x[sample$order[seq_len(h)], , drop = FALSE]
            center <- colMeans(initial)
            scatter <- stats::var(initial) * constants$c1[j] * constants$c2[j]
            inverse <- solve(scatter)
            centered <- sweep(x, 2L, center, "-")
            md <- rowSums((centered %*% inverse) * centered)
            if (j > 1L) {
                error <- abs(md - (h - 1) * quadratic / scale) / pmax(1, abs(md))
                refresh_error <- max(refresh_error, error)
            }
            precision <- inverse * scale / (h - 1)
            quadratic <- md * scale / (h - 1)
        } else {
            md <- (h - 1) * quadratic / scale
        }
        if (any(!is.finite(md)) || any(md < -1e-8))
            stop("Invalid squared distance in recursive update.")
        distance[j] <- sort(md, partial = h + 1L)[h + 1L]
        if (j < length(sizes)) {
            # Update the inverse scatter sum and all centered quadratic forms.
            delta <- x[sample$order[h + 1L], ] - center
            projection <- drop(precision %*% delta)
            norm2 <- sum(delta * projection)
            offset <- drop(x %*% projection) - sum(center * projection)
            weight <- h / (h + 1)
            denominator <- 1 + weight * norm2
            if (!is.finite(denominator) || denominator <= 0)
                stop("Invalid denominator in the inverse scatter update.")
            quadratic <- quadratic - 2 * offset / (h + 1) + norm2 / (h + 1)^2 -
                weight / denominator * (offset - norm2 / (h + 1))^2
            precision <- precision - (weight / denominator) * tcrossprod(projection)
            center <- center + delta / (h + 1)
        }
    }
    path <- data.frame(h = sizes, d_next_sq = distance,
                       threshold_sq = constants$threshold_sq)
    choice <- window_choice(distance, sample$radius)
    fs <- signal_index(path)
    before <- seq_len(fs)
    near_crossing <- abs(distance[before] - path$threshold_sq[before]) <=
        1e-6 * pmax(1, abs(distance[before]))
    list(path = path, ah = choice$index, fs = fs,
         use_direct = !is.finite(refresh_error) || refresh_error > 1e-7 ||
             choice$near_tie || any(near_crossing))
}

adaptive_fits <- function(sample, both, cache) {
    result <- tryCatch(recursive_path(sample, cache), error = function(e) NULL)
    at_indices <- function(path, ah, fs) {
        indices <- if (both) c(ah, fs) else ah
        lapply(path$h[indices], function(h) fit_subset(sample, h, path))
    }
    if (!is.null(result) && !result$use_direct) {
        fits <- tryCatch(at_indices(result$path, result$ah, result$fs),
                         error = function(e) NULL)
        if (!is.null(fits)) {
            exact <- vapply(fits, function(f) rank_distance(f$raw_distance, f$h),
                            numeric(1))
            indices <- if (both) c(result$ah, result$fs) else result$ah
            approximate <- result$path$d_next_sq[indices]
            errors <- abs(exact - approximate) / pmax(1, abs(exact))
            if (all(is.finite(errors)) && all(errors <= 1e-7)) return(fits)
        }
    }
    path <- direct_path(sample)
    at_indices(path, window_choice(path$d_next_sq, sample$radius)$index,
               signal_index(path))
}

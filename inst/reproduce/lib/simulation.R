simulation_design <- function() {
  types <- c("point", "random", "cluster", "point_mass")
  clean <- data.frame(study = "clean", setting = c("A", "B", "C"),
                      type = "clean", eps = 0, r = 0)
  fixed <- expand.grid(type = types, setting = c("A", "B", "C"),
                       eps = c(0.1, 0.3), stringsAsFactors = FALSE)
  fixed$study <- "fixed"
  fixed$r <- 15
  curve <- expand.grid(r = seq(6, 10, 0.5), type = types,
                       eps = c(0.1, 0.3), stringsAsFactors = FALSE)
  curve$study <- "curve"
  curve$setting <- "B"
  design <- rbind(clean, fixed[names(clean)], curve[names(clean)])
  design$id <- seq_len(nrow(design))
  design$n <- unname(c(A = 100L, B = 500L, C = 2000L)[design$setting])
  design$p <- unname(c(A = 10L, B = 50L, C = 200L)[design$setting])
  design
}

latent_sample <- function(config) {
  n <- config$n
  p <- config$p
  zero <- numeric(p)
  identity <- diag(p)
  y <- MASS::mvrnorm(n, zero, identity)
  if (config$type == "clean") return(y)

  if (config$type %in% c("point", "point_mass")) {
    w <- stats::rchisq(n, df = 3)
    y <- sqrt(1 / 3) * sweep(y, 1L, sqrt(w / 3), "/")
  }
  m <- round(n * config$eps)
  indices <- sample.int(n, size = m)
  r <- config$r
  if (config$type == "point") {
    v <- as.numeric(sample(0:(10 * p), p, replace = TRUE))
    v <- v - sum(v) / p
    u <- v / sqrt(sum(v^2))
    y[indices, ] <- MASS::mvrnorm(m, r * sqrt(p) * u, 0.1^2 * identity)
  } else if (config$type == "point_mass") {
    v <- stats::rnorm(p)
    u <- v / sqrt(sum(v^2))
    y[indices, ] <- matrix(r * sqrt(p) * u, m, p, byrow = TRUE)
  } else if (config$type == "random") {
    radius <- r * p^(1 / 4)
    for (i in indices) {
      v <- MASS::mvrnorm(1, zero, identity)
      center <- (radius / sqrt(sum(v^2))) * v
      y[i, ] <- MASS::mvrnorm(1, center, identity)
    }
  } else if (config$type == "cluster") {
    center <- rep(r * p^(1 / 4) / sqrt(p), p)
    y[indices, ] <- MASS::mvrnorm(m, center, identity)
  } else {
    stop("Unknown contamination design: ", config$type)
  }
  y
}

losses <- function(fit, inverse_g) {
  p <- nrow(inverse_g)
  center <- as.vector(inverse_g %*% fit$center)
  scatter <- inverse_g %*% fit$cov %*% inverse_g
  eigenvalues <- Mod(eigen(scatter, symmetric = FALSE, only.values = TRUE)$values)
  c(e_mu = sqrt(sum(center^2)),
    e_Sigma = log10(max(eigenvalues) / min(eigenvalues)),
    MSE_Sigma = sum((scatter - diag(p))^2) / p^2,
    L_Stein = sum(diag(scatter)) - log(det(scatter)) - p)
}

simulate_batch <- function(job) {
  cfg <- job$config
  g <- matrix(0.75, cfg$p, cfg$p)
  diag(g) <- 1
  inverse_g <- if (cfg$study != "curve") solve(g)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  rows <- lapply(job$rep, function(i) {
    set.seed(620000000L + cfg$id * 10000L + i)
    x <- latent_sample(cfg) %*% g
    depth_seed <- 720000000L + cfg$id * 10000L + i
    fits <- FDBAH::fdb_four(x, seed = depth_seed)
    if (cfg$study == "curve") {
      return(data.frame(id = cfg$id, rep = i, h = fits[["FDB-AH"]]$h))
    }
    if (cfg$study == "fixed") fits <- fits[c("FDB-AH", "FDB-FS", "FDB-75")]
    values <- t(vapply(fits, losses, numeric(4L), inverse_g = inverse_g))
    colnames(values) <- metric_names
    if (any(!is.finite(values))) stop("Nonfinite loss: design ", cfg$id, ", replicate ", i)
    data.frame(id = cfg$id, rep = i, method = names(fits), values, row.names = NULL)
  })
  do.call(rbind, rows)
}

simulate <- function(design, map) {
  blocks <- split(1:1000, rep(1:40, each = 25))
  jobs <- lapply(seq_len(nrow(design)), function(j) {
    lapply(blocks, function(ids) list(config = design[j, , drop = FALSE], rep = ids))
  })
  do.call(rbind, map(do.call(c, jobs), simulate_batch))
}

metric_names <- c("e_mu", "e_Sigma", "MSE_Sigma", "L_Stein")
type_labels <- c(point = "Point", point_mass = "Point mass", random = "Random", cluster = "Cluster")

paper_number <- function(x) {
  if (x != 0 && (abs(x) < 0.001 || abs(x) >= 1000)) {
    return(sub("e([+-])0+", "e\\1", sprintf("%.3e", x)))
  }
  sprintf("%.3f", x)
}

table_entries <- function(raw, id, method, metric) {
  rows <- raw[raw$id == id & raw$method == method, c("rep", metric)]
  stopifnot(identical(sort(rows$rep), 1:1000))
  x <- rows[[metric]]
  c(mean = paper_number(mean(x)), sd = paper_number(stats::sd(x)))
}

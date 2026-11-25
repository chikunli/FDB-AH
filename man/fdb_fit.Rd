\name{fdb_fit}
\alias{fdb_fit}
\alias{fdb_four}
\title{Select an Initial Subset and Estimate Multivariate Location and Scatter}
\description{
Rank observations by random projection depth, choose an initial subset, and
apply one step of reweighting. FDB-AH selects the first window maximum of
the ordered squared Mahalanobis distances.
}
\usage{
fdb_fit(x, method = "FDB-AH", seed = 1L, directions = NULL,
        window = 0.08)

fdb_four(x, seed = 1L, timed = FALSE, directions = NULL,
         window = 0.08)
}
\arguments{
  \item{x}{Numeric matrix or data frame. Rows are observations and columns
  are variables. Missing rows are omitted; at least \eqn{p+2} complete rows
  are required for \eqn{p} variables.}
  \item{method}{One of \code{"FDB-AH"}, \code{"FDB-FS"}, \code{"FDB-HF"},
  and \code{"FDB-75"}.}
  \item{seed}{Positive integer seed passed to the random projection algorithm.}
  \item{directions}{Number of random directions. The default is
  \code{max(1000, 10 * ncol(x))}.}
  \item{window}{Window radius as a fraction of the sample size.
  The integer radius is \code{max(1, floor(window * n))}.}
  \item{timed}{If \code{FALSE}, use the paper's shared-path simulation
  implementation. If \code{TRUE}, use the real-data implementation and
  measure elapsed seconds per method, including the common depth computation.}
}
\details{
For a depth-ordered subset of size \eqn{h}, the path value is the
\eqn{(h+1)}st smallest squared distance among all observations. The
covariance estimate uses the Gaussian MCD consistency and finite-sample
corrections. The recursion refreshes the inverse scatter every 64 steps.
Near ties, unstable updates, or discrepancies at a selected size trigger
direct recomputation from subset covariances.

FDB-FS uses the first crossing of the Gaussian reference envelope;
FDB-HF and FDB-75 use fixed sizes. The nominal squared-distance cutoff
for reweighting is the 0.975 quantile of the chi-squared distribution
with \eqn{p} degrees of freedom. In the shared-path simulation the
fixed-size controls use a strict comparison; the real-data routines
include the cutoff. A fit with at most \eqn{p} retained observations
returns its initial location and scatter estimates.
}
\value{
\code{fdb_fit} returns a list with the initial size \code{h},
the original row indices \code{initial} and \code{retained}, the reweighted
\code{center} and \code{cov}, and the initial \code{raw_center} and
\code{raw_cov}. The vector \code{raw_distance} contains squared distances
in complete-row order. For an adaptive method, \code{path} is a data frame
with candidate sizes \code{h}, order statistics \code{d_next_sq}, and
thresholds \code{threshold_sq}; it is \code{NULL} for a fixed-size method.

\code{fdb_four} returns a named list of four fits using the same projection
ordering. With \code{timed = TRUE}, each fit also has an \code{elapsed} field,
and FDB-FS stops evaluating its path at the first crossing.
}

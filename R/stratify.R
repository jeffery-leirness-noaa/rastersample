#' Stratify
#'
#' Stratify a vector or `SpatRaster` object into the given number of strata.
#'
#' @param x Either a numeric vector or a `SpatRaster` object. If `SpatRaster`,
#' only the first layer will be used.
#' @param n_strata Number of strata.
#' @param equal_split Whether or not to stratify into equally-spaced quantiles
#' of the data values. If `TRUE` (the default), `probs` and `vals` should both
#' be `NULL`.
#' @param probs A numeric vector of probabilities, sorted increasingly, for the
#' cut points with values in [0,1]. The length of this vector should equal
#' `n_strata` - 1.
#' @param vals A numeric vector of break point values, sorted increasingly. The
#' length of this vector should equal `n_strata` - 1.
#'
#' @return An object of the same type as `x` (either a numeric vector or a
#' `SpatRaster` object).
#' @export
stratify <- function(x, n_strata, equal_split = TRUE, probs = NULL, vals = NULL) {
  if (equal_split & (!is.null(probs) | !is.null(vals))) {
    stop("If `equal_split`, both `probs` and `vals` must be NULL.")
  }
  if (!equal_split & is.null(probs) & is.null(vals)) {
    stop("If NOT `equal_split`, either `probs` or `vals` must be specified (i.e., must be non-NULL).")
  }
  if (!equal_split & !is.null(probs) & !is.null(vals)) {
    stop("Both `probs` and `vals` are non-NULL.")
  }
  if (!equal_split & !is.null(probs) & length(probs) != (n_strata - 1)) {
    stop("Length of `probs` should equal the number of splits (i.e., `n_strata` - 1).")
  }
  if (!equal_split & !is.null(vals) & length(vals) != (n_strata - 1)) {
    stop("Length of `vals` should equal the number of splits (i.e., `n_strata` - 1).")
  }
  if (inherits(x, what = "SpatRaster")) {
    x <- terra::subset(x, subset = 1)
    df <- x |>
      terra::as.data.frame(cells = TRUE, na.rm = TRUE) |>
      tibble::as_tibble()
    vec <- dplyr::pull(df, var = -1)
  }
  if (equal_split) {
    breaks <- stats::quantile(vec, probs = seq(0, 1, length.out = n_strata + 1))
  } else {
    if (!is.null(probs)) {
      breaks <- stats::quantile(vec, probs = c(0, probs, 1))
    } else if (!is.null(vals)) {
      breaks <- c(min(vec), vals, max(vec))
    }
  }
  if (inherits(x, what = "SpatRaster")) {
    y <- findInterval(vec, vec = breaks, rightmost.closed = TRUE)
    terra::set.values(x, cells = df$cell, values = y)
    names(x) <- paste0(names(x), "_strata")
    x
  } else {
    findInterval(x, vec = breaks, rightmost.closed = TRUE)
  }
}

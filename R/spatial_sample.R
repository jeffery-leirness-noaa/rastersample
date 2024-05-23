#' Spatial sample
#'
#' Take a sample of the data (without replacement) according to the specified
#' sampling method
#'
#' @param x A `data.frame` or `SpatRaster` object.
#' @param n The number of `data.frame` rows or `SpatRaster` cells to select
#' (i.e., the sample size).
#' @param method Which sampling method should be used to select the rows?
#'   * `"random"` (the default): random sampling via `dplyr::slice_sample()`.
#'   * `"biased"`: random sampling of rows for which `bias_var` value is greater
#'   than `bias_thresh` via `dplyr::slice_sample()`.
#'   * `"stratified"`: stratified random sampling. Randomly select (n / number
#'   of strata) rows for each value of `strata_var` via `dplyr::slice_sample()`.
#'   * `"clh"`: conditioned latin hypercube sampling via `clhs::clhs()`.
#'   * `"balanced"`: spatially balanced sampling via `MBHdesign::quasiSamp()`.
#'   * `"balanced-stratified"`: spatially balanced stratified sampling via
#'   `MBHdesign::quasiSamp()`.
#' @param bias_var The `data.frame` column or `SpatRaster` layer to use for
#' biased sampling (i.e., when `method` is `"biased"`). Must be a character
#' vector of length 1.
#' @param bias_thresh The biased sampling numeric threshold value. If `method`
#' is `"biased"`, only rows for which `bias_var` value is greater than
#' `bias_thresh` are considered eligible for sampling.
#' @param clh_var The `data.frame` columns or `SpatRaster` layers to use for
#' conditioned latin hypercube sampling (i.e., when `method` is `"clh"`). Must
#' be character vector of length greater than 1.
#' @param clh_iter A positive number, giving the number of iterations for the
#' Metropolis-Hastings annealing process. If `method` is `"clh"`, this value is
#' passed to the `iter` argument of `clhs::clhs()`.
#' @param strata_var The `data.frame` column(s) or `SpatRaster` layer(s) that
#' define the strata for stratified sampling (i.e., when `method` is
#' `"stratafied"` or `"balanced-stratified"`).
#' @param drop_na Whether or not to exclude `data.frame` rows or `SpatRaster`
#' cells with NA values prior to sampling.
#' @param as_raster Whether or not to return a `SpatRaster` object.
#'
#' @return Either a `data.frame` (tibble) object (if `as_raster` is `FALSE`) or
#' a `SpatRaster` object (if `as_raster` is `TRUE`).
#' @export
spatial_sample <- function(x, n, method, bias_var = NULL, bias_thresh = NULL, clh_var = NULL, clh_iter = NULL, strata_var = NULL, drop_na = TRUE, as_raster = FALSE) {

  if (method %in% c("balanced", "balanced-stratified")) {
    stopifnot("if `method` is 'balanced' or 'balanced-stratified', then x must be a SpatRaster" = inherits(x, "SpatRaster"))
    if (method == "balanced-stratified") {
      stopifnot("if `method` is 'balanced-stratified', then `strata_var` must be specified" = !is.null(strata_var))
    }
  }

  if (drop_na) {
    if (inherits(x, "SpatRaster")) {
      x[!terra::noNA(x)] <- NA
    } else {
      df <- x |>
        tidyr::drop_na() |>
        tibble::as_tibble()
    }
  }

  if (inherits(x, "SpatRaster") & !(method %in% c("balanced", "balanced-stratified"))) {
    df <- x |>
      terra::as.data.frame(cells = TRUE, na.rm = drop_na) |>
      tibble::as_tibble()
  }

  if (method == "random") {
    samp <- spatialsample_random(df, n = n)
  } else if (method == "biased") {
    samp <- spatialsample_biased(df, n = n, var = bias_var, thresh = bias_thresh)
  } else if (method == "stratified") {
    samp <- spatialsample_stratified(df, n = n, var = strata_var)
  } else if (method == "clh") {
    samp <- spatialsample_clh(df, n = n, var = clh_var, iter = clh_iter)
  } else if (method %in% c("balanced", "balanced-stratified")) {
    samp <- spatialsample_balanced(x, n = n, strata_var = strata_var, drop_na = drop_na)
  }

  if (inherits(x, "SpatRaster") & as_raster) {
    x_samp <- x
    x_samp[-samp$cell] <- NA
    x_samp
  } else {
    samp
  }

}
# Spatial sample of data
spatialsample_random <- function(data, n) {
  data |>
    tibble::as_tibble() |>
    dplyr::slice_sample(n = n)
}
spatialsample_biased <- function(data, n, var, thresh) {
  data |>
    tibble::as_tibble() |>
    dplyr::filter(.data[[var]] > thresh) |>
    dplyr::slice_sample(n = n)
}
spatialsample_stratified <- function(data, n, var) {
  n_strata <- data |>
    dplyr::select(dplyr::all_of(var)) |>
    dplyr::n_distinct()
  data |>
    tibble::as_tibble() |>
    dplyr::group_by(.data[[var]]) |>
    dplyr::slice_sample(n = ceiling(n / n_strata))
}
spatialsample_clh <- function(data, n, var, iter) {
  clh_samp <- data |>
    dplyr::select(dplyr::all_of(var)) |>
    clhs::clhs(size = n, iter = iter, use.cpp = TRUE, simple = FALSE,
               progress = TRUE)
  data |>
    tibble::as_tibble() |>
    dplyr::slice(clh_samp$index_samples)
}
spatialsample_balanced <- function(r, n, strata_var = NULL, drop_na = TRUE) {
  if (!is.null(strata_var)) {
    df <- r |>
      terra::as.data.frame(cells = TRUE, na.rm = drop_na) |>
      tibble::as_tibble() |>
      dplyr::group_by(.data[[strata_var]]) |>
      dplyr::mutate(.w = 1 / dplyr::n()) |>
      dplyr::ungroup()
    r_ip <- r |>
      terra::subset(subset = strata_var)  # converting NA's to 0's (as in `balanced` sample above) causes an issue with the raster resulting from terra::set.values() below
    terra::set.values(r_ip, cells = df$cell, values = df$.w * 1e+07)
  } else {
    if (drop_na) {
      r_ip <- terra::noNA(r) * 1
    } else {
      r_ip <- terra::subset(r, subset = 1)
      terra::values(r_ip) <- 1
    }
  }
  ba_samp <- MBHdesign::quasiSamp(n = n, dimension = 2, inclusion.probs = r_ip) |>
    tibble::as_tibble()
  r |>
    terra::as.data.frame(cells = TRUE, na.rm = drop_na) |>
    tibble::as_tibble() |>
    dplyr::filter(cell %in% ba_samp$ID)
}

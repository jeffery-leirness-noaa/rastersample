#' Spatial sample
#'
#' Take a sample of the data according to the specified sampling design
#'
#' @param data A `data.frame` object, ideally the output from
#' spatialsample_prep().
#' @param n The number of `data.frame` rows to select (i.e., the sample size).
#' @param design Which sampling design should be used to select the rows?
#'   * `"random"` (the default): random sampling via `dplyr::slice_sample()`.
#'   * `"biased"`: random sampling of rows for which `bias_var` value is greater than
#'   `bias_thresh` via `dplyr::slice_sample()`.
#'   * `"stratified"`: stratified random sampling. Randomly select (n / number of strata) rows for each
#'   value of `strata_var` via `dplyr::slice_sample()`.
#'   * `"clh"`: conditioned latin hypercube sampling via `clhs::clhs()`.
#'   * `"balanced"`: spatially balanced sampling via `MBHdesign::quasiSamp()`.
#'   * `"balanced_stratified"`: spatially balanced stratified sampling via
#'   `MBHdesign::quasiSamp()`.
#' @param r A `SpatRaster` object specifying the
#' @param bias_var description
#' @param bias_thresh description
#' @param strata_var description
#' @param clh_var description
#' @param clh_iter A positive number, giving the number of iterations for the
#' Metropolis-Hastings annealing process. If `design` is `"clh"`, this value is
#' passed to the `iter` argument of `clhs::clhs()`.
#'
#' @return description
#' @export
spatial_sample <- function(data, n, design, r, bias_var, bias_thresh, strata_var, clh_var, clh_iter) {
  if (design == "random") {
    spatialsample_random(data, n = n)
  } else if (design == "biased") {
    spatialsample_biased(data, n = n, var = bias_var, thresh = bias_thresh)
  } else if (design == "stratified") {
    spatialsample_stratified(data, n = n, var = strata_var)
  } else if (design == "clh") {
    spatialsample_clh(data, n = n, var = clh_var, iter = clh_iter)
  } else if (design == "balanced") {
    spatialsample_balanced(data, n = n, r = r)
  } else if (design == "balanced_stratified") {
    spatialsample_balanced_stratified(data, n = n, r = r, var = strata_var)
  }
}
#' Spatial sample of data
#'
#' @inheritParams spatial_sample
#' @param var description
#'
#' @return description
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
spatialsample_balanced <- function(data, n, r) {
  r_ip <- r |>
    terra::subset(subset = 1) |>
    terra::not.na() * 1
  ba_samp <- MBHdesign::quasiSamp(n = n, dimension = 2, inclusion.probs = r_ip) |>
    tibble::as_tibble()
  data |>
    tibble::as_tibble() |>
    dplyr::filter(.cell %in% ba_samp$ID)
}
spatialsample_balanced_stratified <- function(data, n, r, var) {
  df <- data |>
    dplyr::group_by(.data[[var]]) |>
    dplyr::mutate(.w = 1 / dplyr::n()) |>
    dplyr::ungroup()
  r_ip <- r |>
    terra::subset(subset = 1)  # converting NA's to 0's (as in `balanced` sample above) causes an issue with the raster resulting from terra::set.values() below
  terra::set.values(r_ip, cells = df$.cell, values = df$.w * 1e+07)
  ba_samp <- MBHdesign::quasiSamp(n = n, dimension = 2, inclusion.probs = r_ip)
  data |>
    tibble::as_tibble() |>
    dplyr::filter(.cell %in% ba_samp$ID)
}

#' Spatial sample prep
#'
#' @param data description
#' @param n description
#' @param design description
#'
#' @return description
spatialsample_prep <- function(x, cells = FALSE, xy = FALSE, scale = FALSE) {
  if (scale) x <- terra::scale(x, center = TRUE, scale = TRUE)
  df <- x |>
    terra::as.data.frame(cells = cells, xy = xy, na.rm = TRUE) |>
    tibble::as_tibble()
  if (cells) {
    df <- df |>
      dplyr::rename(.cell = cell) |>
      dplyr::relocate(.cell, .after = dplyr::last_col())
  }
  if (xy) {
    df <- df |>
      dplyr::rename(.x = x, .y = y) |>
      dplyr::relocate(.x, .y, .after = dplyr::last_col())
  }
  df
}

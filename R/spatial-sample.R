#' Spatial sample
#'
#' Take a sample of the data (without replacement) according to the specified
#' sampling method. This function provides multiple sampling strategies for spatial data,
#' ranging from simple random sampling to more sophisticated spatially balanced approaches.
#'
#' @param x A `data.frame` or `SpatRaster` object. For spatial sampling methods
#'   ('balanced', 'balanced-stratified') and line sampling, a `SpatRaster` is required.
#' @param n The number of `data.frame` rows or `SpatRaster` cells to select
#'   (i.e., the sample size).
#' @param method Which sampling method should be used to select the rows?
#'   * `"random"` (default): random sampling via `dplyr::slice_sample()`.
#'   * `"biased"`: random sampling of rows for which `bias_var` value is greater
#'   than `bias_thresh` via `dplyr::slice_sample()`.
#'   * `"stratified"`: stratified random sampling. Randomly select (n / number
#'   of strata) rows for each value of `strata_var` via `dplyr::slice_sample()`.
#'   * `"clh"`: conditioned latin hypercube sampling via `clhs::clhs()`.
#'   * `"balanced"`: spatially balanced sampling via `MBHdesign::quasiSamp()`.
#'   * `"balanced-stratified"`: spatially balanced stratified sampling via
#'   `MBHdesign::quasiSamp()`.
#' @param bias_var The `data.frame` column or `SpatRaster` layer to use for
#'   biased sampling (i.e., when `method` is `"biased"`). Must be a character
#'   vector of length 1.
#' @param bias_thresh The biased sampling numeric threshold value. If `method`
#'   is `"biased"`, only rows for which `bias_var` value is greater than
#'   `bias_thresh` are considered eligible for sampling.
#' @param clh_var The `data.frame` columns or `SpatRaster` layers to use for
#'   conditioned latin hypercube sampling (i.e., when `method` is `"clh"`). Must
#'   be character vector of length greater than 1.
#' @param clh_iter A positive number, giving the number of iterations for the
#'   Metropolis-Hastings annealing process. If `method` is `"clh"`, this value is
#'   passed to the `iter` argument of `clhs::clhs()`.
#' @param strata_var The `data.frame` column(s) or `SpatRaster` layer(s) that
#'   define the strata for stratified sampling (i.e., when `method` is
#'   `"stratified"` or `"balanced-stratified"`).
#' @param drop_na Whether or not to exclude `data.frame` rows or `SpatRaster`
#'   cells with NA values prior to sampling. Default is `TRUE`.
#' @param as_raster Whether or not to return a `SpatRaster` object. Default is `FALSE`.
#' @param type The type of sampling to perform. Options are:
#'   * `"point"` (default): sample individual points/cells/rows.
#'   * `"line"`: sample transect lines (requires `x` to be a `SpatRaster` and
#'     `method` to be either `"random"` or `"balanced"`).
#' @param control A list of additional parameters to pass to
#'   `MBHdesign::transectSamp()`. This argument is only used when `type` is
#'   `"line"`.
#'
#' @return Either a `data.frame` (tibble) object (if `as_raster` is `FALSE`) or
#'   a `SpatRaster` object (if `as_raster` is `TRUE`).
#'   - For point sampling, the returned data frame contains the sampled points with
#'     their associated values and cell IDs if from a raster.
#'   - For line sampling, the returned data frame contains the sampled transect
#'     points with their coordinates and cell values.
#'
#' @examples
#' \dontrun{
#' library(terra)
#'
#' # Create a sample raster
#' r <- rast(nrows = 100, ncols = 100)
#' values(r) <- 1:ncell(r)
#'
#' # Simple random sampling
#' samples <- spatial_sample(r, n = 10, method = "random")
#'
#' # Spatially balanced sampling
#' balanced_samples <- spatial_sample(r, n = 10, method = "balanced")
#'
#' # Stratified sampling using a categorical raster
#' r_strata <- r
#' values(r_strata) <- sample(LETTERS[1:5], ncell(r), replace = TRUE)
#' strat_samples <- spatial_sample(r, n = 50, method = "stratified",
#'                                 strata_var = "lyr.1")
#'
#' # Line sampling
#' transects <- spatial_sample(r, n = 5, method = "balanced", type = "line")
#' }
#'
#' @seealso
#' \code{\link[dplyr]{slice_sample}} for the underlying random sampling function.
#' \code{\link[clhs]{clhs}} for conditioned latin hypercube sampling.
#' \code{\link[MBHdesign]{quasiSamp}} for spatially balanced sampling.
#' \code{\link[MBHdesign]{transectSamp}} for transect sampling.
#'
#' @importFrom dplyr filter slice_sample group_by mutate ungroup n_distinct select
#' @importFrom terra as.data.frame subset noNA extract values set.values
#' @importFrom tibble as_tibble
#' @importFrom tidyr drop_na
#' @importFrom clhs clhs
#' @importFrom MBHdesign quasiSamp transectSamp
#'
#' @export
spatial_sample <- function(
  x,
  n,
  method,
  bias_var = NULL,
  bias_thresh = NULL,
  clh_var = NULL,
  clh_iter = NULL,
  strata_var = NULL,
  drop_na = TRUE,
  as_raster = FALSE,
  type = "point",
  control = NULL
) {
  if (method %in% c("balanced", "balanced-stratified")) {
    stopifnot(
      "if `method` is 'balanced' or 'balanced-stratified', then x must be a SpatRaster" = inherits(
        x,
        "SpatRaster"
      )
    )
    if (method == "balanced-stratified") {
      stopifnot(
        "if `method` is 'balanced-stratified', then `strata_var` must be specified" = !is.null(
          strata_var
        )
      )
    }
  }
  if (type == "line") {
    stopifnot(
      "if `type` is 'line', then x must be a SpatRaster" = inherits(
        x,
        "SpatRaster"
      )
    )
    stopifnot(
      "if `type` is 'line', then `method` must be either 'random' or 'balanced'" = method %in%
        c("random", "balanced")
    )
  }

  if (drop_na) {
    if (inherits(x, "SpatRaster")) {
      x[!terra::noNA(x)] <- NA
    } else {
      df <- x |>
        tidyr::drop_na()
    }
  } else if (!inherits(x, "SpatRaster")) {
    if (inherits(x, "sf")) {
      df <- x
    } else {
      df <- tibble::as_tibble(x)
    }
  }

  if (
    inherits(x, "SpatRaster") &
      !(method %in% c("balanced", "balanced-stratified"))
  ) {
    df <- x |>
      terra::as.data.frame(cells = TRUE, na.rm = drop_na) |>
      tibble::as_tibble()
  }

  if (type == "line") {
    samp <- spatialsample_transect(
      x,
      n = n,
      method = method,
      drop_na = drop_na,
      control = control
    )
  } else {
    if (method == "random") {
      samp <- spatialsample_random(df, n = n)
    } else if (method == "biased") {
      samp <- spatialsample_biased(
        df,
        n = n,
        var = bias_var,
        thresh = bias_thresh
      )
    } else if (method == "stratified") {
      samp <- spatialsample_stratified(df, n = n, var = strata_var)
    } else if (method == "clh") {
      samp <- spatialsample_clh(df, n = n, var = clh_var, iter = clh_iter)
    } else if (method %in% c("balanced", "balanced-stratified")) {
      samp <- spatialsample_balanced(
        x,
        n = n,
        strata_var = strata_var,
        drop_na = drop_na
      )
    }
  }

  if (inherits(x, "SpatRaster") & as_raster) {
    x_samp <- x
    x_samp[-samp$cell] <- NA
    x_samp
  } else {
    samp
  }
}

#' @rdname spatial_sample
#' @keywords internal
#' @description
#' `spatialsample_random()` performs simple random sampling from a data frame.
spatialsample_random <- function(data, n) {
  data |>
    dplyr::slice_sample(n = n)
}

#' @rdname spatial_sample
#' @keywords internal
#' @description
#' `spatialsample_biased()` performs biased sampling where only rows meeting a
#' threshold criterion are considered.
spatialsample_biased <- function(data, n, var, thresh) {
  data |>
    dplyr::filter(.data[[var]] > thresh) |>
    dplyr::slice_sample(n = n)
}

#' @rdname spatial_sample
#' @keywords internal
#' @description
#' `spatialsample_stratified()` performs stratified random sampling based on a stratification variable.
spatialsample_stratified <- function(data, n, var) {
  n_strata <- data |>
    tibble::as_tibble() |>
    dplyr::select(dplyr::all_of(var)) |>
    dplyr::n_distinct()
  data |>
    dplyr::group_by(.data[[var]]) |>
    dplyr::slice_sample(n = ceiling(n / n_strata))
}

#' @rdname spatial_sample
#' @keywords internal
#' @description
#' `spatialsample_clh()` performs conditioned latin hypercube sampling.
spatialsample_clh <- function(data, n, var, iter) {
  if (is.null(iter)) {
    iter <- 1e+04
  }
  clh_samp <- data |>
    dplyr::select(dplyr::all_of(var)) |>
    clhs::clhs(
      size = n,
      iter = iter,
      use.cpp = TRUE,
      simple = FALSE,
      progress = TRUE
    )
  data |>
    dplyr::slice(clh_samp$index_samples)
}

#' @rdname spatial_sample
#' @keywords internal
#' @description
#' `spatialsample_balanced()` performs spatially balanced sampling on a SpatRaster.
#' If `strata_var` is provided, it performs stratified spatially balanced sampling.
spatialsample_balanced <- function(r, n, strata_var = NULL, drop_na = TRUE) {
  if (!is.null(strata_var)) {
    df <- r |>
      terra::as.data.frame(cells = TRUE, na.rm = drop_na) |>
      tibble::as_tibble() |>
      dplyr::group_by(.data[[strata_var]]) |>
      dplyr::mutate(.w = 1 / dplyr::n()) |>
      dplyr::ungroup()
    r_ip <- r |>
      terra::subset(subset = strata_var) # converting NA's to 0's (as in `balanced` sample above) causes an issue with the raster resulting from terra::set.values() below
    terra::set.values(r_ip, cells = df$cell, values = df$.w * 1e+07)
  } else {
    if (drop_na) {
      r_ip <- terra::noNA(r) * 1
    } else {
      r_ip <- terra::subset(r, subset = 1)
      terra::values(r_ip) <- 1
    }
  }
  ba_samp <- MBHdesign::quasiSamp(
    n = n,
    dimension = 2,
    inclusion.probs = r_ip
  ) |>
    tibble::as_tibble()
  r |>
    terra::as.data.frame(cells = TRUE, na.rm = drop_na) |>
    tibble::as_tibble() |>
    dplyr::filter(cell %in% ba_samp$ID)
}

#' @rdname spatial_sample
#' @keywords internal
#' @description
#' `spatialsample_transect()` performs transect-based sampling on a SpatRaster.
spatialsample_transect <- function(
  r,
  n,
  method,
  drop_na = TRUE,
  control = NULL
) {
  if (drop_na) {
    r_ip <- terra::subset(r, subset = 1)
    r_ip[!is.na(r_ip)] <- 1
  } else {
    r_ip <- terra::subset(r, subset = 1)
    terra::values(r_ip) <- 1
  }
  names(r_ip) <- "values"
  df <- terra::as.data.frame(r_ip, xy = TRUE, na.rm = FALSE) |>
    tibble::as_tibble() |>
    dplyr::arrange(y, x)

  spat_random_type <- dplyr::case_match(
    method,
    "random" ~ "pseudo",
    "balanced" ~ "quasi"
  )
  control_transect <- list(
    transect.pattern = "line",
    spat.random.type = spat_random_type
  )
  if (!is.null(control)) {
    control_transect <- append(control_transect, control)
  }

  ba_samp <- MBHdesign::transectSamp(
    n,
    potential.sites = as.matrix(dplyr::select(df, x, y)),
    inclusion.probs = df$values,
    control = control_transect
  )
  cells <- terra::extract(
    r,
    dplyr::select(ba_samp$points, x, y),
    cells = TRUE,
    ID = FALSE
  )
  dplyr::bind_cols(ba_samp$points, cells) |>
    tibble::as_tibble()
}

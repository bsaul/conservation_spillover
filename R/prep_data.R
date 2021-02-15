
#' Computes the length of shared boundaries between a parcel and its neighbors
#' 
#' @param x the `MULTIPOLYGON` for the parcel of interest
#' @param neighbors a list of `MULTIPOLYGON`s of which share a boundary with `x`
#' @keywords  internal
.shared_boundaries_length <- function(x, neighbors){
  purrr::map(
    .x = neighbors,
    .f = ~ {
      sf::st_intersection(x, .x) %>%
        sf::st_cast(warn = FALSE) %>%
        sf::st_length() %>%
        sum()
    }
  )
}

#' Computes the length of shared boundaries between a parcel and its neighbors
#' 
#' @param df a spatial dataframe, containing a variable named by `boundary_col` 
#'          (defaulting to `"boundary"`) which is a `LINES` object
#' @param i the row index in `df` of the parcel of interest
#' @param nghrs the indices of neighbors of `df[i, ]` in `df`
#' @returns a `units` vector containing boundary lengths in meters
#' @export
shared_boundaries_length <- function(df, i, ngbrs, boundary_col = "boundary"){
  `if`(
    length(ngbrs) > 0,
    .shared_boundaries_length(
      x = df[i, ], 
      neighbors = df[ngbrs, boundary_col])[[boundary_col]],
    units::set_units(numeric(), "m")
  )
}
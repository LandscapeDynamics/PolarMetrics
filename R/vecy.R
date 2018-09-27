#' @name vec.y
#' @title Calculate vertical component vectors
#' @description Calculate the vertical vector component for each pair
#'   of polar coordinates.
#' @param r A vector of numeric values representing dates in radians.
#' @param v A vector of numeric values representing a uniformly sampled
#'   series of data that vary over multiple cycles.
#' @details vec.y calculates the length of the vertical vector component
#'   of each polar coordinate pair. Given an angle \eqn{r} in radians
#'   and a corresponding amplitude \eqn{v} the vertical vector component,
#'   \eqn{Vy}, is calculated by:
#' \deqn{VY(r, v) = v * \sin(r)}{VY(r, v) = v * sin(r)}.
#' @return Returns the length of the vertical vector component for each
#'   polar coordinate pair.
#' @examples
#' dpy <- 365                 # Days/yr
#' data(mndvi)                # Load data
#' t <- as.vector(mndvi$day)  # Days since January 1, 2000
#' r <- t2rad(t,dpy)          # Transform days of year to radians
#' v <- as.vector(mndvi$wc)   # MODIS NDVI for Willow Creek tower, WI
#' vy <- vec.y(r,v)           # The vertical vector components
#' @author Bjorn J. Brooks, Danny C. Lee, William W. Hargrove, Lars Y. Pomara
#' @references Brooks, B.J., Lee, D.C., Desai, A.R., Pomara, L.Y.,
#'   Hargrove, W.W. (2017). Quantifying seasonal patterns in
#'   disparate environmental variables using the PolarMetrics R package.
#' @export

vec.y <- function(r,v) {
  r <- as.numeric(r)                               # Ensure type is numeric
  v <- as.numeric(v)
  if (length(r) == length(v)) {
    return(v * sin(r))                             # Vert vec comp. of r,v
  } else {
    stop('Number of values in arg 1 should = arg 2')
  }
}

#' @name t2rad
#' @title Transform dates to radians
#' @description Transform dates to radians.
#' @param t A vector of numeric dates to be transformed to polar coordinates.
#' @param dpc A numeric value representing the number of divisions per
#'   cycle (e.g., hours per day, days per year, etc.)
#' @details \code{t2rad} transforms the sequence \eqn{t} of dates to
#'   radians using:
#' \deqn{\frac{(t - 1) \bmod d}{d} * 2\pi}{(t - 1) \% d / d * 2\pi}
#' where \eqn{d} indicates the number of days per year in sequence \eqn{t}.
#' @return Returns a vector of polar coordinates the same length as
#'   input vector \eqn{t}.
#' @examples
#' dpy <- 365                 # Days/year
#' data(mndvi)                # Load data
#' t <- as.vector(mndvi$day)  # Days since January 1, 2000
#' r <- t2rad(t, dpc=dpy)     # Transform days of year to radians
#' head(cbind(t,r))           # Compare results
#' @author Bjorn J. Brooks, Danny C. Lee, William W. Hargrove, Lars Y. Pomara
#' @references Brooks, B.J., Lee, D.C., Desai, A.R., Pomara, L.Y.,
#'   Hargrove, W.W. (2017). Quantifying seasonal patterns in
#'   disparate environmental variables using the PolarMetrics R package.
#' @export

t2rad <- function(t, dpc) {
  if (length(dpc) == 1 & length(t) >= length(dpc)) {
    return(2 * pi * (t %% dpc) / dpc)              # Series t in radians
  } else {
    stop('Number of values in arg 1 should be > arg 2')
  }
}

#' @name rad2d
#' @title Transform from radians to time units
#' @description Transform radians to dates in time units (e.g., day of year).
#' @param r A vector of numeric values of angles in radians (representing
#'   times/dates) to be converted back to time units.
#' @param dpc A numeric value representing the number of divisions per
#'   cycle (e.g., hours per day, days per year, etc.)
#' @details \code{rad2d} transforms input values of \code{r} from radians
#'   to dates in the interval [0,dpc+1).
#' @return Returns a vector of cyclic dates in the interval [0,dpc+1)
#'   having the same length as input vector \eqn{r}.
#' @examples
#' dpy <- 365                 # Days/year
#' data(mndvi)                # Load data
#' t <- as.vector(mndvi$day)  # Days since January 1, 2000
#' r <- t2rad(t, dpc=dpy)     # Transform days of year to radians
#' days_from_rads <- rad2d(r, dpc=dpy) # Transform radians back to days of yr
#' head(cbind(t, days_from_rads)) # Compare results. Cols should be equiv.
#' @author Bjorn J. Brooks, Danny C. Lee, William W. Hargrove, Lars Y. Pomara
#' @references Brooks, B.J., Lee, D.C., Desai, A.R., Pomara, L.Y.,
#'   Hargrove, W.W. (2017). Quantifying seasonal patterns in
#'   disparate environmental variables using the PolarMetrics R package.
#' @export

rad2d <- function(r, dpc) {
  if (length(dpc) == 1 & length(r) >= length(dpc)) {
    return(as.integer(dpc * r / (2 * pi)))         # Transf. units to days
  } else {
    stop('arg 1 should have 1 value and arg 2 >= 1 value')
  }
}

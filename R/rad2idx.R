#' @name rad2idx
#' @title Transform radians to sequential index values
#' @description Transform radian values to the nearest corresponding indices of the original data.
#' @param r A vector of numeric values representing dates in radians.
#' @param spc A numeric value indicating the number of samples per cycle (e.g., samples per year).
#' @details rad2idx finds the corresponding indices for each element of \eqn{r}.
#' @return Returns the index values corresponding to \eqn{r}.
#' @examples
#' dpy <- 365                 # Days/year
#' spy <- 46                  # Samples/year (or samples/cycle)
#' data(mndvi)                # Load data
#' t <- as.vector(mndvi$day)  # Days since January 1, 2000
#' r <- t2rad(t, dpc=dpy)     # Transform days of year to radians
#' r_idx <- rad2idx(r, spc=spy) # Indices corresponding to radial angles in r
#' head(cbind(t,r,r_idx))     # Compare results
#' @author Bjorn J. Brooks, Danny C. Lee, William W. Hargrove, Lars Y. Pomara
#' @references Brooks, B.J., Lee, D.C., Desai, A.R., Pomara, L.Y.,
#'   Hargrove, W.W. (accepted). Quantifying seasonal patterns in
#'   disparate environmental variables using the PolarMetrics R package.
#' @export

rad2idx <- function(r, spc) {
  if (length(spc) == 1 & length(r) >= length(spc)) {
    c <- spc/(2*pi) # Coefficient for conversion from radians to index
    idx1 <- r*c # Convert to index units
    idx2 <- ceiling(idx1) # Round up to nearest whole index value

    return(idx2)
  } else {
    stop('Arg 2 should have 1 element & arg 1 >= 1 element')
  }
}

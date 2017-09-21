#' @name vec_mag
#' @title Calculate Euclidean magnitudes from horizontal and vertical
#'   component vectors
#' @description Calculate the magnitude (Euclidean distance) for each
#'   pair of horizontal and vertical component vectors.
#' @param vx A numeric vector of 1 or more values representing
#'   horizontal components of Euclidean vectors. \code{vx} can be obtained
#'   using \code{\link{vec.x}}.
#' @param vy A numeric vector of 1 or more values representing vertical
#'  components of Euclidean vectors. \code{vy} can be obtained using
#'  \code{\link{vec.y}}.
#' @details \code{vec_mag} returns the Euclidean distance of each pair
#'   of horizontal and vertical component vectors according to:
#' \deqn{\sqrt{vx^2 + vy^2}}{\sqrt{vx^2 + vy^2}}
#' @return Returns the Euclidean distance corresponding to each horizontal
#'   and vertical component vector pair.
#' @examples
#' dpy <- 365                 # Days/yr
#' data(mndvi)                # Load data
#' t <- as.vector(mndvi$day)  # Days since January 1, 2000
#' r <- t2rad(t,dpy)          # Transform days of year to radians
#' v <- as.vector(mndvi$wc)   # MODIS NDVI for Willow Creek tower, WI
#' vx <- mean(vec.x(r,v), na.rm=TRUE) # Avg horizontal vector
#' vy <- mean(vec.y(r,v), na.rm=TRUE) # Avg vertical vector
#' vm <- vec_mag(vx,vy)       # Magnitude (length) of average vector
#' @author Bjorn J. Brooks, Danny C. Lee, William W. Hargrove, Lars Y. Pomara
#' @references Brooks, B.J., Lee, D.C., Desai, A.R., Pomara, L.Y.,
#'   Hargrove, W.W. (accepted). Quantifying seasonal patterns in
#'   disparate environmental variables using the PolarMetrics R package.
#' @export

vec_mag <- function(vx,vy) {
  if (length(vx) == length(vy)) {
    return(sqrt(vx^2+vy^2))	# Length of vector
  } else {
    stop('Number of values in arg 1 should = arg 2')
  }
}

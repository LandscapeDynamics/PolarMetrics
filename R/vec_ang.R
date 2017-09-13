#' @name vec_ang
#' @title Calculate radial angles from horizontal and vertical component vectors
#' @description Calculate the radial angle for each horizontal and vertical component vector pair.
#' @param vx A numeric vector of 1 or more values representing horizontal components of Euclidean vectors. \code{vx} can be obtained using \code{\link{vec.x}}.
#' @param vy A numeric vector of 1 or more values representing vertical components of Euclidean vectors. \code{vy} can be obtained using \code{\link{vec.y}}.
#' @details vec_ang returns the angle corresponding to the horizontal (\eqn{vx}) and vertical vector components (\eqn{vy}) according to:
#' \deqn{atan2(vy, vx) + \pi}{atan2(vy, vx) + \pi}
#' @return Returns the polar coordinate angle corresponding to each horizontal and vertical component vector pair.
#' @examples
#' dpy <- 365                 # Days/yr
#' data(mndvi)                # Load data
#' t <- as.vector(mndvi$day)  # Days since January 1, 2000
#' r <- t2rad(t,dpy)          # Transform days of year to radians
#' v <- as.vector(mndvi$wc)   # MODIS NDVI for Willow Creek tower, WI
#' vx <- mean(vec.x(r,v), na.rm=TRUE) # Avg horizontal vector
#' vy <- mean(vec.y(r,v), na.rm=TRUE) # Avg vertical vector
#' rv_ang <- vec_ang(vx,vy)    # Angle of resultant vec (point of max activity)
#' @author Bjorn J. Brooks, Danny C. Lee, William W. Hargrove, Lars Y. Pomara
#' @references Brooks, B.J., Lee, D.C., Desai, A.R., Pomara, L.Y.,
#'   Hargrove, W.W. (accepted). Quantifying seasonal patterns in
#'   disparate environmental variables using the PolarMetrics R package.
#' @export

vec_ang <- function(vx,vy) {
  if (length(vx) == length(vy)) {
    return(atan2(vy,vx)%%(2*pi))	# vector of angles (0->2 pi rads)
  } else {
    stop('Number of values in arg 1 should = arg 2')
  }
}

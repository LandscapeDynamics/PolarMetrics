#' @name avec_ang
#' @title Calculate the opposites of the input radial angle(s)
#' @description Calculate the vertical opposite of the input radian value(s).
#' @param r A numeric vector of values of radial coordinates.
#' @details avec_ang returns the vertical opposite(s) of the input angles.
#' @return Returns the vertical opposite of radial angle(s).
#' @examples
#' dpy <- 365                 # Days/yr
#' data(mndvi)                # Load data
#' t <- as.vector(mndvi$day)  # Days since January 1, 2000
#' r <- t2rad(t,dpy)          # Transform days of year to radians
#' v <- as.vector(mndvi$wc)   # MODIS NDVI for Willow Creek tower, WI
#' vx <- mean(vec.x(r,v), na.rm=TRUE) # Avg horizontal vector
#' vy <- mean(vec.y(r,v), na.rm=TRUE) # Avg vertical vector
#' vm <- vec_mag(vx,vy)       # Magnitude of average vector
#' rv_ang <- vec_ang(vx,vy)   # Angle of resultant vec (point of max activity)
#' av_ang <- avec_ang(rv_ang)  # Angle marking point of least activity
#' @author Bjorn J. Brooks, Danny C. Lee, William W. Hargrove, Lars Y. Pomara
#' @references Brooks, B.J., Lee, D.C., Desai, A.R., Pomara, L.Y.,
#'   Hargrove, W.W. (2017). Quantifying seasonal patterns in
#'   disparate environmental variables using the PolarMetrics R package.
#' @export

###
avec_ang <- function(r) {
  return((r + pi) %% (2 * pi))	                    # Opp. of r (0->2 pi rads)
}

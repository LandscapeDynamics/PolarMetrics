#' @name sum_cycle
#' @title Calculate cumulative sums within each cycle in a series
#' @description Calculate the cumulative sums within each cycle in a uniformly
#'   sampled series. This function is useful where the phenological (or
#'   variable-defined) year is offset from the calendar year. Note that 1) This
#'   function assumes every cycle contains the same number of elements, and 2)
#'   This function truncates the values before the start of the first cycle and
#'   after the end of the last cycle, which for a series of \eqn{n} cycles
#'   results in an output series with \eqn{n - 1} cycles.
#' @param v A vector of numeric values representing a uniformly sampled series
#'   of data that oscillate over multiple cycles.
#' @param b A positive integer indicating the index (i.e., 1 -> \eqn{spc}) in
#'   the average cycle corresponding to the beginning (or period of least
#'   activity) of the average cycle. Within polar space, \code{b} lies along a
#'   ray that intersects the origin and evenly bisects the summed values of
#'   average cycle. \code{b} can be obtained using \code{\link{avec_ang}}.
#' @param spc A numeric value indicating the number of samples per cycle (e.g., samples per year).
#' @details \code{sum_cycle} calculates the cumulative sums within cycles
#'   beginning from a specified offset (i.e., argument \eqn{b > 1}). Let \eqn{b}
#'   be the index for the beginning of cycle \eqn{c} and \eqn{e} be the index
#'   for the ending index for the cycle \eqn{\{v_{b}, v_{b+1}, \dots,
#'   v_{e}\}}{\{v_b, v_b+1, \dots, v_e\}}. Within each cycle the sum at each
#'   index \eqn{j} is calculated as the sum of all values leading up to and
#'   including \eqn{j} following: \deqn{c_j = \sum_{i=b}^j v_i, v_{i+1}, \dots
#'   v_{j}}{c_j = v_i + v_i+1 + \dots + v_j} This is repeated over \eqn{c - 1}
#'   cycles in series \eqn{v} to obtain a new series of cumulative sums.
#' @return Returns a vector of values (and their corresponding indices from the
#'   input vector) representing the cumulative sums within each cycle based on a
#'   specified beginning point. If \code{b} > 1 then the output vector will have
#'   \code{spc} fewer values, which were truncated due to the offset. Output
#'   vector will begin at \code{b}, the average minimum not at the beginning of
#'   the time series, unless the average min is the beginning (\code{b} = 1).
#' @examples
#' dpy <- 365                 # Days/yr
#' spy <- 46                  # Number of samples in one cycle (yr)
#' data(mndvi)                # Load data
#' t <- as.vector(mndvi$day)  # Days since January 1, 2000
#' r <- t2rad(t, dpc=dpy)     # Transform days of year to radians
#' v <- as.vector(mndvi$wc)   # MODIS NDVI for Willow Creek tower, WI
#' vx <- mean(vec.x(r,v), na.rm=TRUE) # Avg horizontal vector
#' vy <- mean(vec.y(r,v), na.rm=TRUE) # Avg vertical vector
#' rv_ang <- vec_ang(vx,vy)  # Angle of the resultant vector (the direction
#'   # that the average vector points)
#' av_ang <- avec_ang(rv_ang)  # Angle marking point of least activity
#' av_idx <- rad2idx(av_ang, spc=spy) # Index (1-spc) marking avg start of yr
#' ann_cum <- sum_cycle(v,av_idx,spc=spy)$cumsum # Accum. vals within each yr
#' @author Bjorn J. Brooks, Danny C. Lee, William W. Hargrove, Lars Y. Pomara
#' @references Brooks, B.J., Lee, D.C., Desai, A.R., Pomara, L.Y.,
#'   Hargrove, W.W. (2017). Quantifying seasonal patterns in
#'   disparate environmental variables using the PolarMetrics R package.
#' @export

sum_cycle <- function(v,b,spc) {
  if (b <= spc & length(spc) == 1 & length(v) %% spc == 0 & length(v) > length(spc) & (length(v)/spc) > 1) {
    c <- length(v)/spc				   # Num. of cycles
    if (b %% spc == 0) { # Special case if natural (phenological) yr not offset
      cs <- c(NA,spc*c)                            # Initialize
      idx <- c(NA,spc*c)                           # Initialize
      for (I in 1:c) {
        begidx.v <- 1+(spc*(I-1))                  # Beginning index
        endidx.v <- spc*I	                   # Ending index
        begidx.cs <- begidx.v
        endidx.cs <- endidx.v
        cs[begidx.cs:endidx.cs] <- cumsum(v[begidx.v:endidx.v]) # Cum sum
        idx[begidx.cs:endidx.cs] <- begidx.v:endidx.v # Corresp. idx in input
      }
    } else { # If natural (phenological) yr is offset
      cs <- c(NA,spc*(c-1)) # Initialize cumsum variable with room for c-1 cycles
      idx <- c(NA,spc*(c-1)) # Initialize index variable with room for c-1 cycles
      for (I in 1:(c-1)) {
        begidx.v <- b+(spc*(I-1))                  # Beginning index
        endidx.v <- b+(spc*I)-1	                   # Ending index
        begidx.cs <- 1+(spc*(I-1))
        endidx.cs <- spc*I
        cs[begidx.cs:endidx.cs] <- cumsum(v[begidx.v:endidx.v]) # Cum sum
        idx[begidx.cs:endidx.cs] <- begidx.v:endidx.v # Corresp. idx in input
      }
    }
  } else if (b > spc) {
    stop('the value of arg 2 should be <= arg 3. See documentation')
  } else if ((length(v)/spc) <= 1) {
    stop('arg 1 vector should contain > 1 cycle/yr of data. See documentation')
  } else if (length(v) %% spc != 0) {
    stop('arg 1 should be exactly divisible by arg 2. See documentation')
  } else {
    stop('arg 1 should have >= 1 value & args 2,3 should have 1 value each')
  }
  output <- data.frame(vidx=idx, cumsum=cs)

  return(output)
}

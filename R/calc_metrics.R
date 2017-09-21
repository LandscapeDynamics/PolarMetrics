#' @name calc_metrics
#' @title Easy function for calculating all polar meatrics. 
#' @description A basic, but fast function for calculating all polar
#'   metrics from arrays or xts objects.
#' @importFrom stats sd
#' @importFrom zoo index
#' @param input Either a vector of values or a column of an xts object
#'   indexed with time stamps. If \code{input} is an xts object then omit
#'   argument \code{t}. Else if \code{input} is a vector then a
#'   corresponding vector of values must be passed to \code{t}.
#' @param t An optional vector of time values (e.g., days) corresponding
#'   to the input vector sampling points. Do not use this argument if
#'   \code{input} is an xts object.
#' @param yr_type Argument specifying either 'cal_yr' for output (of
#'   timing variables) given in days starting from Jan. 1, or 'offset_yr'
#'   for output in days starting from the average seasonal minimum.
#' @param spc Integer value specifying the number of samples per cycle
#'   (measurements per year) in input.
#' @param lcut Numeric value in the range [0,0.5] passed to
#'   \code{window_idx} function. Indicates the percentile to truncate
#'   from the first half of each cycle. For example, 0.15 will result
#'   in remove the interval corresponding to [0\%,15\%] representing a
#'   window that begins after the 15th-percentile is crossed.
#' @param hcut Numeric value in the range (0.5,1] passed to
#'   \code{window_idx} function. Indicates the percentile to truncate
#'   from the last half of each cycle. For example, 0.85 will result in
#'   remove the interval corresponding to [85\%,100\%] representing a
#'   window that begins after the 85th-percentile is crossed.
#' @param return.vecs logical argument specifying whether or not to
#'   include all of the horizontal and vertical component vectors in output.
#' @details \code{calc_metrics} runs through the entire polar
#'   transformation process and conveniently outputs the final polar
#'   metrics for all years included in the input.
#' @return Returns a list with the derived polar metrics (e.g., length of
#'   season). Optional objects in the list are the actual vecturs used
#'   in deriving the polar metrics and the average vector angle and
#'   magnitude.
#' @examples
#' library(PolarMetrics)
#' library(xts)
#' input <- xts(mndvi$fef, as.Date(mndvi$date))
#' calc_metrics(input, yr_type='cal_yr', spc=46, lcut=0.15, hcut=0.8,
#'              return.vecs=FALSE) # Calculate metrics
#' ### Return the average vectors for the entire time series 
#' calc_metrics(input, yr_type='cal_yr', spc=46, lcut=0.15, hcut=0.8,
#'              return.vecs=TRUE)$avg.vectors
#' ### Return the horizontal and vertical vector components
#' head(calc_metrics(input, yr_type='cal_yr', spc=46, lcut=0.15, hcut=0.8,
#'              return.vecs=TRUE)$vectors)
#' @author Bjorn J. Brooks, Danny C. Lee, William W. Hargrove, Lars Y. Pomara
#' @references Brooks, B.J., Lee, D.C., Desai, A.R., Pomara, L.Y.,
#'   Hargrove, W.W. (accepted). Quantifying seasonal patterns in
#'   disparate environmental variables using the PolarMetrics R package.
#' @export

calc_metrics <- function(input, t=NULL, yr_type, spc, lcut, hcut, return.vecs) {
  dpy <- 365 # Days per year
  if (is.null(t)) { # If no values supplied for t then obtain from xts index
    t <- as.numeric(format(index(input), '%j')) # Time coordinates
  }
  t <- as.vector(t)
  input <- as.vector(input)
  if (length(input) != length(t)) {
    stop('input and t should be have the same number of values')
  }
  nyr <- length(input)/spc  # No. of years in input
  r <- t2rad(t, dpc=dpy)       # Transform days of year to radians
  v <- as.vector(input)
  VX <- vec.x(r,v)        # Horizontal vectors
  VY <- vec.y(r,v)        # Vertical vectors
  vx <- mean(VX, na.rm=T) # Avg horizontal vector
  vy <- mean(VY, na.rm=T) # Avg vertical vector
  rv_ang <- vec_ang(vx,vy) # Angle of resultant vector (point of max activity)
  rv_idx <- rad2idx(rv_ang, spc=spc) # Index (1-spc) marking avg vector
  rv_doy <- rad2d(rv_ang, dpc=dpy) # Day of year equivalent of rv_ang
  rv_mag <- vec_mag(vx,vy) # Magnitude (length) of resultant vector
  avec_ang <- avec_ang(rv_ang) # Angle marking point of least activity
  avec_doy <- rad2d(avec_ang, dpc=dpy) # Day of year equivalent of avec_ang
  avec_idx <- rad2idx(avec_ang, spc=spc) # Index (1-spc) marking avg start of yr
  ann_cum <- sum_cycle(v,avec_idx,spc=spc)$cumsum # Accum vals within offset yrs
  npy <- length(ann_cum)/spc # No. of complete years in output
  # Note, output is re-centered, and accumulates starting from avec_idx
  #  (offset of variable yr from calendar yr) and has nyr-1 yrs of data
  #  due to re-centering.

  output <- data.frame(yr=rep(NA,npy), es=rep(NA,npy), ems=rep(NA,npy),
                    ms=rep(NA,npy), lms=rep(NA,npy), ls=rep(NA,npy),
                    s_intv=rep(NA,npy), s_avg=rep(NA,npy),
                          s_sd=rep(NA,npy), s_mag=rep(NA,npy),
                          ems_mag=rep(NA,npy), lms_mag=rep(NA,npy))
  for (J in 1:npy) { # Calculate pheno param's for each yr
    output$yr[J] <- J
    wi <- window_idx(ann_cum,npy,J,lcut,hcut) # early, mid, late-S ann_cum indx
    es_idx <- wi[1]
    ems_idx <- wi[2]
    ms_idx <- wi[3]
    lms_idx <- wi[4]
    ls_idx <- wi[5]
    es <- t[es_idx]   # Day of yr marking lcut (e.g. 15 %tile) threshold
    ems <- t[ems_idx] # Day of yr marking lcut+50/2 threshold
    ms <- t[ms_idx]   # Day of yr marking 50 %tile threshold
    lms <- t[lms_idx] # Day of yr marking hcut+50/2 threshold
    ls <- t[ls_idx]   # Day of yr marking hcut threshold
    if (yr_type == 'cal_yr') {
      # Timing variables relative to the calendar year
      output$es[J] <- round((es + avec_doy) %% dpy)   # DOY for ES milestone
      output$ems[J] <- round((ems + avec_doy) %% dpy) # DOY for EMS milestone
      output$ms[J] <- round((ms + avec_doy) %% dpy)   # DOY for MS milestone
      output$lms[J] <- round((lms + avec_doy) %% dpy) # DOY for MLS milestone
      output$ls[J] <- round((ls + avec_doy) %% dpy)   # DOY for LS milestone
    } else if (yr_type == 'offset_yr') {
      # Timing variables relative to the offset (variable-centered) year
      output$es[J] <- es %% dpy   # DOY marking ES milestone
      output$ems[J] <- ems %% dpy # DOY marking EMS milestone
      output$ms[J] <- ms %% dpy   # DOY marking MS milestone
      output$lms[J] <- lms %% dpy # DOY marking MLS milestone
      output$ls[J] <- ls %% dpy   # DOY marking LS milestone
    }
    output$s_intv[J] <- ls - es # Days in the growing season
    output$s_avg[J] <- mean(v[es_idx:ls_idx], na.rm=TRUE) # Mean for Seasn.
    output$s_sd[J] <- sd(v[es_idx:ls_idx]) # Standard deviation for Seasn.
    # Magnitude (length) of average vector
    output$s_mag[J] <- vec_mag(mean(VX[es_idx:ls_idx], na.rm=TRUE),
                              mean(VY[es_idx:ls_idx], na.rm=TRUE))
    # Magnitude of avg vec between ES & MS thresholds
    output$ems_mag[J] <- vec_mag(mean(VX[es_idx:(ms_idx-1)],
                                      na.rm=TRUE),
                          mean(VY[es_idx:(ms_idx-1)], na.rm=TRUE))
    # Magnitude of avg vec between MS & LS thresholds
    output$lms_mag[J] <- vec_mag(mean(VX[ms_idx:(ls_idx-1)],
                                      na.rm=TRUE),
                          mean(VY[ms_idx:(ls_idx-1)], na.rm=TRUE))
  }

  if (return.vecs == FALSE) {
    return(output)
  }  else if (return.vecs == TRUE) {
    # If vector components requested then return by combining as list object
    output <- list('metrics' = output, 'vectors' = data.frame(VX=VX, VY=VY),
                   'avg.vectors' = data.frame(rv_idx=rv_idx, rv_ang=rv_ang,
                                   rv_doy=rv_doy, rv_mag=rv_mag,
                                   avec_idx=avec_idx, avec_ang=avec_ang,
                                   avec_doy=avec_doy))
    return(output)
  }
}

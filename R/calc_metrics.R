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
#'   timing variables) given in days starting from Jan. 1, or 'rot_yr'
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
#' @param return_vecs logical argument specifying whether or not to
#'   include all of the horizontal and vertical component vectors in output.
#' @param sin_cos logical argument. If TRUE then each timing
#'   metric (es, ms, etc.) is returned as its sine and cosine components,
#'   es is returned as es_sin and es_cos.
#' @details \code{calc_metrics} runs through the entire polar
#'   transformation process and conveniently outputs the final polar
#'   metrics for all years included in the input.
#' @return Returns a list with all of the derived polar metrics (e.g.,
#'   early season, mid season, etc.). Timing variables in output can
#'   be returned relative to the standard calendar year or rotated to
#'   a relative year using \code{yr_type} argument.
#'   Sine and cosine components can also be returned instead of days since
#'   start using \code{sin_cos} argument.
#'   The optional argument \code{return_vecs} can be used to add
#'   the actual vecturs used in deriving the polar metrics to the returned
#'   list. This will also return the overall average (all years) statistics
#'   of the resultant vector (rv) and its opposite the anti-vector (av).
#'   Table below indicates the variable components of each object within
#'   the returned list.
#' metrics (data frame):     year,
#'                           es (or es_sin, es_cos),    # Early ssn DOY
#'                           ems (or ems_sin, ems_cos), # Early-mid ssn DOY
#'                           ms (or ms_sin, ms_cos),    # Mid ssn DOY
#'                           lms (or lms_sin, lms_cos), # Late-mid ssn DOY
#'                           ls (or ls_sin, ls_cos),    # Late ssn DOY
#'                           s_intv,                    # Season length(days)
#'                           s_avg,                     # Avg data val in ssn
#'                           s_sd,                      # StDev of data dur ssn
#'                           s_mag,                     # Mag of avg vec in ssn
#'                           ems_mag,                   # Mag early-mid ssn vec
#'                           lms_mag,                   # Mag late-mid ssn vec
#'                           a_avg                      # Avg data val of yr
#' vectors (data frame):     VX, VY
#' avg_vectors (data frame): rv_idx, rv_ang, rv_doy, rv_mag,
#'                           avec_idx, avec_ang, avec_doy
#' @examples
#' library(PolarMetrics)
#' library(xts)
#' input <- xts(mndvi$fef, as.Date(mndvi$date))
#' ### Calculate polar measures relative to calendar year
#' calc_metrics(input, yr_type='cal_yr', spc=46, lcut=0.15, hcut=0.8,
#'              return_vecs=FALSE, sin_cos=FALSE)
#' ### Calculate as above and return sine, cosine components of timing metrics
#' calc_metrics(input, yr_type='cal_yr', spc=46, lcut=0.15, hcut=0.8,
#'              return_vecs=FALSE, sin_cos=TRUE)
#' ### Calculate & return the average vectors for the entire time series
#' calc_metrics(input, yr_type='cal_yr', spc=46, lcut=0.15, hcut=0.8,
#'              return_vecs=TRUE, sin_cos=FALSE)$avg_vectors
#' ### Calculate & return the horizontal and vertical vector components
#' head(calc_metrics(input, yr_type='cal_yr', spc=46, lcut=0.15, hcut=0.8,
#'              return_vecs=TRUE, sin_cos=FALSE)$vectors)
#' @author Bjorn J. Brooks, Danny C. Lee, William W. Hargrove, Lars Y. Pomara
#' @references Brooks, B.J., Lee, D.C., Desai, A.R., Pomara, L.Y.,
#'   Hargrove, W.W. (accepted). Quantifying seasonal patterns in
#'   disparate environmental variables using the PolarMetrics R package.
#' @export

calc_metrics <- function(input, t=NULL, yr_type, spc, lcut, hcut, return_vecs, sin_cos) {
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
  ann_cum <- sum_cycle(v,avec_idx,spc=spc)$cumsum # Accum vals within rot. yrs
  npy <- length(ann_cum)/spc # No. of complete years in output
  # Note, output is re-centered, and accumulates starting from avec_idx
  #  (offset of variable yr from calendar yr) and has nyr-1 yrs of data
  #  due to re-centering.

  if (!isTRUE(sin_cos)) { # if false
    output <- data.frame(yr=rep(NA,npy), es=rep(NA,npy), ems=rep(NA,npy),
                         ms=rep(NA,npy), lms=rep(NA,npy), ls=rep(NA,npy),
                         s_intv=rep(NA,npy), s_avg=rep(NA,npy),
                         s_sd=rep(NA,npy), s_mag=rep(NA,npy),
                         ems_mag=rep(NA,npy), lms_mag=rep(NA,npy),
                         a_avg=rep(NA,npy))
  } else if (isTRUE(sin_cos)) { # if true
    output <- data.frame(yr=rep(NA,npy),
			 es_sin=rep(NA,npy), es_cos=rep(NA,npy),
			 ems_sin=rep(NA,npy), ems_cos=rep(NA,npy),
			 ms_sin=rep(NA,npy), ms_cos=rep(NA,npy),
			 lms_sin=rep(NA,npy), lms_cos=rep(NA,npy),
			 ls_sin=rep(NA,npy), ls_cos=rep(NA,npy),
                         s_intv=rep(NA,npy), s_avg=rep(NA,npy),
                         s_sd=rep(NA,npy), s_mag=rep(NA,npy),
                         ems_mag=rep(NA,npy), lms_mag=rep(NA,npy),
			                   a_avg=rep(NA,npy))
  }

  for (J in 1:npy) { # Calculate pheno param's for each yr
    output$yr[J] <- J
    wi <- window_idx(ann_cum,npy,J,lcut,hcut) # early, mid, late-S ann_cum indx
    es_idx <- wi[1]
    ems_idx <- wi[2]
    ms_idx <- wi[3]
    lms_idx <- wi[4]
    ls_idx <- wi[5]
    be_idx = c((spc*(J-1)+1),(spc*J)) # This cycle's beg/end indices
    es <- t[es_idx]   # Day of yr marking lcut (e.g. 15 %tile) threshold
    ems <- t[ems_idx] # Day of yr marking lcut+50/2 threshold
    ms <- t[ms_idx]   # Day of yr marking 50 %tile threshold
    lms <- t[lms_idx] # Day of yr marking hcut+50/2 threshold
    ls <- t[ls_idx]   # Day of yr marking hcut threshold
    if (yr_type == 'cal_yr' & !isTRUE(sin_cos)) {
      # Timing variables relative to the calendar year
      output$es[J] <- round((es + avec_doy) %% dpy)   # DOY for ES milestone
      output$ems[J] <- round((ems + avec_doy) %% dpy) # DOY for EMS milestone
      output$ms[J] <- round((ms + avec_doy) %% dpy)   # DOY for MS milestone
      output$lms[J] <- round((lms + avec_doy) %% dpy) # DOY for MLS milestone
      output$ls[J] <- round((ls + avec_doy) %% dpy)   # DOY for LS milestone
    } else if (yr_type == 'rot_yr' & !isTRUE(sin_cos)) {
      # Timing variables relative to the rotated year
      output$es[J] <- es %% dpy   # DOY marking ES milestone
      output$ems[J] <- ems %% dpy # DOY marking EMS milestone
      output$ms[J] <- ms %% dpy   # DOY marking MS milestone
      output$lms[J] <- lms %% dpy # DOY marking MLS milestone
      output$ls[J] <- ls %% dpy   # DOY marking LS milestone
    } else if (yr_type == 'cal_yr' & isTRUE(sin_cos)) {
      # Sine, cosine timing variables relative to the calendar year
      output$es_sin[J] <- sin((es + avec_doy) %% dpy)   # sin(DOY) for ES
      output$es_cos[J] <- cos((es + avec_doy) %% dpy)   # cos(DOY) for ES
      output$ems_sin[J] <- sin((ems + avec_doy) %% dpy) # sin(DOY) for EMS
      output$ems_cos[J] <- cos((ems + avec_doy) %% dpy) # cos(DOY) for EMS
      output$ms_sin[J] <- sin((ms + avec_doy) %% dpy)   # sin(DOY) for MS
      output$ms_cos[J] <- cos((ms + avec_doy) %% dpy)   # cos(DOY) for MS
      output$lms_sin[J] <- sin((lms + avec_doy) %% dpy) # sin(DOY) for LMS
      output$lms_cos[J] <- cos((lms + avec_doy) %% dpy) # cos(DOY) for LMS
      output$ls_sin[J] <- sin((ls + avec_doy) %% dpy)   # sin(DOY) for LS
      output$ls_cos[J] <- cos((ls + avec_doy) %% dpy)   # cos(DOY) for LS
    } else if (yr_type == 'rot_yr' & isTRUE(sin_cos)) {
      # Sine, cosine timing vars relative to the rotated yr
      output$es_sin[J] <- sin(es %% dpy)   # sin(DOY) for ES
      output$es_cos[J] <- cos(es %% dpy)   # cos(DOY) for ES
      output$ems_sin[J] <- sin(ems %% dpy) # sin(DOY) for EMS
      output$ems_cos[J] <- cos(ems %% dpy) # cos(DOY) for EMS
      output$ms_sin[J] <- sin(ms %% dpy)   # sin(DOY) for MS
      output$ms_cos[J] <- cos(ms %% dpy)   # cos(DOY) for MS
      output$lms_sin[J] <- sin(lms %% dpy) # sin(DOY) for LMS
      output$lms_cos[J] <- cos(lms %% dpy) # cos(DOY) for LMS
      output$ls_sin[J] <- sin(ls %% dpy)   # sin(DOY) for LS
      output$ls_cos[J] <- cos(ls %% dpy)   # cos(DOY) for LS
    }
    output$s_intv[J] <- ls - es # Days in the growing season
    output$s_avg[J] <- mean(v[es_idx:ls_idx], na.rm=TRUE) # Mean for Seasn.
    output$s_sd[J] <- sd(v[es_idx:ls_idx]) # Standard deviation for Seasn.
    output$a_avg[J] <- mean(v[be_idx[1]:be_idx[2]], na.rm=TRUE) # Full yr mean
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

  if (!isTRUE(return_vecs)) { # if false
    return(output)
  }  else if (isTRUE(return_vecs)) { # if true
    # If vector components requested then return by combining as list object
    output <- list('metrics' = output, 'vectors' = data.frame(VX=VX, VY=VY),
                   'avg_vectors' = data.frame(rv_idx=rv_idx, rv_ang=rv_ang,
                                   rv_doy=rv_doy, rv_mag=rv_mag,
                                   avec_idx=avec_idx, avec_ang=avec_ang,
                                   avec_doy=avec_doy))
    return(output)
  }
}

\encoding{UTF-8}
\name{PolarMetrics-package}
\alias{PolarMetrics}
\docType{package}
\keyword{package}
\title{Functions for analyzing seasonal patterns}
\description{
This package can be used on any uniformly sampled environmental time
series containing a periodic cycle (i.e., temperature, CO2, NDVI,
stream flow). These functions transform and reproject time series on
to a polar plane in order to calculate the natural, i.e. data-defined,
seasonal timing dates, seasonal asymmetry and magnitudes. Timing metrics
can be expressed relative to calendar dates or based on a "natural" year
that begins at the point of least activity (shortest vectors). Although
developed to improve interpretation and analysis for NDVI time series
these functions work on other types of data that have strong periodicities
of any length (daily, annual, interannual).
}

\details{
\code{PolarMetrics} is a family of primarily trigonometric functions
written for extracting annual polar measures from uniformly sampled
environmental time series. These functions ultimately output a set of
timing and magnitude measures that describe attributes of each year
(cycle), such as seasonality (deciduousness). These measures are useful
for subsequent analysis such as clustering and factor analysis where
specific selections of polar metrics (e.g., length of season, seasonal
magnitude) are ingested to explore particular features or saesons. The
\code{window_idx} function excludes user-defined portions of the start and
end of each cycle, which isolates the primary season of interest. In the
output, timing metrics can be given relative to the calendar year (days
1-365 corresponding to Jan. 1-Dec. 31), or the offset year (days 1-365
starting from the average point of least activity). Although uncommon,
note that if the natural year and calendar years are syncrhonized the
polar metrics output will have the same number of years as the input
data. However, in most cases the natural year in the data is offset and
thus the output measures will have one less year of data than in the
input, because the values preceeding the first offset are discarded as
well as the values in the final year following the last complete year.\cr

Table 1. Example output of derived polar metrics for MODIS NDVI from
2000-2010 for a single pixel covering the Coweeta Hydrologic Lab eddy
flux tower in North Carolina. (See list of variable definitions below)
\tabular{llllllllllll}{
  Year \tab  ES \tab EMS \tab MS  \tab LMS \tab LS  \tab Sintv \tab Savg \tab Ssd  \tab Smag  \tab EMSma \tab LMSma\cr
  2000 \tab  98 \tab 162 \tab 210 \tab 258 \tab 306 \tab 208 \tab 0.702 \tab 0.126 \tab 0.397 \tab 0.536 \tab 0.699\cr
  2001 \tab  98 \tab 170 \tab 218 \tab 258 \tab 306 \tab 208 \tab 0.704 \tab 0.120 \tab 0.390 \tab 0.522 \tab 0.730\cr
  2002 \tab 106 \tab 170 \tab 218 \tab 266 \tab 314 \tab 208 \tab 0.726 \tab 0.111 \tab 0.406 \tab 0.571 \tab 0.709\cr
  2003 \tab 106 \tab 162 \tab 218 \tab 258 \tab 306 \tab 200 \tab 0.750 \tab 0.094 \tab 0.431 \tab 0.601 \tab 0.739\cr
  2004 \tab 106 \tab 162 \tab 218 \tab 258 \tab 306 \tab 200 \tab 0.727 \tab 0.116 \tab 0.425 \tab 0.572 \tab 0.727\cr
  2005 \tab 106 \tab 170 \tab 218 \tab 266 \tab 306 \tab 200 \tab 0.727 \tab 0.129 \tab 0.431 \tab 0.562 \tab 0.743\cr
  2006 \tab 106 \tab 170 \tab 218 \tab 258 \tab 314 \tab 208 \tab 0.733 \tab 0.120 \tab 0.415 \tab 0.580 \tab 0.718\cr
  2007 \tab  98 \tab 162 \tab 218 \tab 258 \tab 306 \tab 208 \tab 0.719 \tab 0.116 \tab 0.401 \tab 0.552 \tab 0.725\cr
  2008 \tab 106 \tab 162 \tab 218 \tab 258 \tab 306 \tab 200 \tab 0.699 \tab 0.113 \tab 0.408 \tab 0.549 \tab 0.701\cr
  2009 \tab 106 \tab 162 \tab 218 \tab 258 \tab 306 \tab 200 \tab 0.723 \tab 0.120 \tab 0.428 \tab 0.577 \tab 0.718\cr
  2010 \tab 106 \tab 162 \tab 218 \tab 258 \tab 306 \tab 200 \tab 0.722 \tab 0.126 \tab 0.429 \tab 0.578 \tab 0.716
}


\subsection{Interpreting polar graphs}{
At first glance polar plots can be difficult to inspect, as compared to
more familiar time series representations. However, polar graphs have
several benefits over linear time traces. By wrapping data points radially
around the center, deviations from the typical timing and progression of
the seasonal cycle are easier to notice. For example the broad width of
orbits in the polar plot in Fig. 1(b) indicate that there is significant
difference among years, but most of this orbital width is constrained
to December through February revealing that most of this interannual
change is due to changes in the winter vegetation signature (evergreen
component of the vegetation). Another visual advantage of polar plots
is that changes in the timing of the season become more obvious. For
example, during September the breadth in orbital path in Fig. 1(b) is
constrained, indicating little interannual variability at that point in
the year. However, during October the path is broader (more variable),
revealing that October is the more variable month.

% Example of polar conversion
\if{html}{\figure{ex1.png}{options: width='653px' alt='Figure: ex1.png'}}
\if{latex}{\figure{ex1.pdf}{options: width=5in}}

Figure 1. Time series and polar plot comparison. Data
are from twelve years of MODIS NDVI measurements for
one pixel (aprox. 250 x 250 m resolution) covering the
\href{http://flux.aos.wisc.edu/twiki/bin/view/Main/ChEASData}{Willow
Creek AmeriFlux tower} in northern Wisconsin. (a) shows a
conventional time series. (b) shows the same data when graphed in
polar coordinates. \code{RVec} (Resultant Vector, calculated by
\code{\link{vec_ang}}) marks the angle of the average polar vector,
which points toward mid-August.  \code{AVec} is its opposite (calculated
by \code{\link{avec_ang}}), which points at the period of least activity
during the average cycle and is also the starting point for the average
phenological year. \code{SMag} (calculated by \code{\link{vec_mag}}),
shows the degree of seasonality (length of red line). The length of
\code{SMag} is small because it is the average displacement of all
vectors from the center of the polar graph, half of which have negative
angles. (This figure produced by code in Example 3 below).
}


\subsection{Interpreting polar metrics}{
Polar metrics are annual quantities that measure the timing and asymmetry
of annual cycles. Consider, for example, an evergreen forest, as measured
by NDVI. Its timing metrics will be more or less equally spaced and its
seasonality will be small. On the other hand a landscape dominated by
deciduous trees, because it has larger differences in greenness between
winter and summer, will have a higher degree of seasonality as well
as timing dates that are skewed closer to the middle of the growing
season. This package calculates metrics that describe timing asymmetry
and amplitude, which allows changes to be compared across years.\cr


\subsection{Calculating timing metrics}{
Threshold completion measures are calculated within each cycle (year) by
starting from the point in the polar year corresponding to the shortest
vector lengths (See \code{AVec} Fig. 1(b) and start of phenological year
in Fig. 2). Values, in this case NDVI, are accumulated within each year
(see \code{sum_cycle}) and each data point is assigned a corresponding
cumulative percentage based on its proportional value within that year's
total. User defined timing thresholds (e.g., 15\% and 80\%) are then
assigned (see \code{window_idx}) as the time (day of year) when the data
first crossed a given threshold. In these examples timing metrics are
measured in days since the start of the calendar year or alternatively the
offset/phenological year. (Phenological year days have the advantage of
controlling for latitudinal, elevational, or other geographical bias in
phenology timing.) Take for example the early season timing metric (ES)
on the phenological year and assuming a 15\% threshold. An \code{ES}
value of 65 days means that 65 days passed until at least 15\% of
the year's total NDVI accumulated. The \code{window_idx} function is
written such that timing metrics can occur on integer value days that
were actually sampled, so the 15\% threshold actually corresponds to
the first sample day when 15\% or more of the year's cumulative total
NDVI was achieved. Fig. 2 provides an illustration of how the threshold
completion measures were calculated for one example year.

% Example of threshold completion milestones
\if{html}{\figure{ex2.png}{options: width='599px' alt='Figure: ex2.png'}}
\if{latex}{\figure{ex2.pdf}{options: width=5in}}

Figure 2. Timing metrics. Data correspond to year two (2001) of NDVI
at Willow creek. Starting from the day corresponding to the average
low point in the seasonal cycle, NDVI is cumulatively summed within the
phenological year. The early season metric (ES) marks the first day when
at least 15\% of cumulative NDVI was reached (day 107 of calendar year
or 65 on phenological year. See example code for \code{\link{window_idx}}
function). (This figure produced by code in Example 4 below) \cr


\describe{
\item{\code{AVec} (start/end of phenological yr)}{Number of days (offset)
between Jan. 1 and the point of least activity (shortest vectors in the
polar plot)}
\item{\code{Early season (ES)}}{Number of days to cross first threshold
in cumulative annual NDVI (e.g., 15\%).}
\item{\code{Early-mid season (EMS)}}{Number of days to cross the midpoint
(e.g., 32.5\%) between first threshold (e.g., 15\%) and midpoint (50\%).}
\item{\code{Mid-season (MS)}}{Number of days required to cross 50\%
threshold}
\item{\code{Mid-late season (LMS)}}{Number of days to cross the midpoint
(e.g., 67.5\%) between mid-season (50\%) and midpoint (e.g., 80\%).}
\item{\code{Late season (LS)}}{Number of days to cross late season
threshold (e.g., 80\%)}
\item{\code{Length of season (LOS)}}{Number of days between early and
late season thresholds}
}
}


\subsection{Calculating asymmetry and amplitude metrics}{
These power measures are calculated from the values falling within the
\code{ES} and \code{LS} thresholds. If \code{\link{window_idx}} is set to
search the entire year (e.g., by setting \code{p1}=0 and \code{p2}=1)
then asymmetry and amplitude measures will be based on the entire
phenological year of data. The list of these measures is as follows:

\describe{
\item{\code{Savg}}{The average NDVI value between early and late
thresholds (e.g., 15-80\%)}
\item{\code{Ssd}}{The standard deviation of NDVI values between early
and late thresholds}
\item{\code{Smag} (seasonality/deciduousness)}{The length (in units of
input data) of the average vector between ES and LS timing points}
\item{\code{EMSmag} ("spring" seasonal magnitude)}{The length of the
average vector between ES and MS timing points}
\item{\code{LMSmag} ("fall" seasonal magnitude)}{The length of the
average vector between MS and LS timing points}
}

% Example of four phenological years plotted separately.
\if{html}{\figure{ex3.png}{options: width='653px' alt='Figure: ex3.png'}}
\if{latex}{\figure{ex3.pdf}{options: width=5in}}

Figure 3. Polar plots of timing and asymmetry metrics for four years
at Fraser Experimental Forest in Colorado. Each panels indicates
annual changes in timing and asymmetry metrics. Note that the angle of
\code{SMag} corresponds to the timing of mid season (MS), but the length
of \code{SMag}, seaonality, is actually a measure of the degree to which
data are massed in the direction of the \code{MS} date. These years in
particular represent the early stages of mortality in lodgepole pines at
FEF caused by bark beetles. Early season NDVI (April) appeared to start a
variable decline in 2001, meanwhile summer NDVI remained fairly constant
(perhaps due to compensatory greening from understory vegetation). Note
that this area had been under significant drought stress and we do not
intend to suggest this as a system for diagnosing specific disturbances,
only that this technique is capable of measuring changes caused by
environmental shifts. (This figure produced by code in Example 5 below)
}

}
}


\author{
Bjorn J. Brooks, Danny C. Lee, William W. Hargrove, Lars Y. Pomara

Maintainer: Bjorn J. Brooks <bjorn@geobabble.org>
}
\references{
Brooks, B.J., Lee, D.C., Desai, A.R., Pomara, L.Y., Hargrove,
W.W. (accepted). Quantifying seasonal patterns in disparate environmental
variables using the PolarMetrics R package.
}

\examples{
### Example 1, (Table 1 above) Use easy function to calculate polar vectors
library(PolarMetrics)
data(mndvi)                # Load data
lb <- 0.15                 # Cutoff point, as a proportion of 1 (0-0.5), to
                           # remove from each years cumulative total.
                           # e.g., 0.1 will remove from 0 to 0.1.
ub <- 0.8                  # Cutoff point, as a proportion of 1 (0.5-1), to
                           # remove at the end of each years cumulative total.
                           # e.g., 0.9 will remove 0.9 to 1.
s_p_yr <- 46               # Number of MODIS values per year (8-day NDVI)
# Calculate polar metrics for one MODIS NDVI pixel within Coweeta Hydrologic Lab
pmets <- calc_metrics(mndvi$chl, t=mndvi$day, yr_type='cal_yr',
  spc=s_p_yr, lcut=lb, hcut=ub, return.vecs=FALSE)

### Example 2, Long version of Example 1
library(PolarMetrics)
pmets <- data.frame(yr=NA, es=NA, ems=NA, ms=NA, lms=NA, ls=NA,
                         s_intv=NA, s_avg=NA, s_sd=NA, s_mag=NA,
                         ems_mag=NA, lms_mag=NA) # Initialize output
dpy <- 365                 # Days/year
lb <- 0.15                 # Lower cutoff for cumulative annual NDVI
ub <- 0.8                  # Upper cutoff for cumulative annual NDVI
c <- 12                    # Num. of years/cycles
cxs <- 1                   # Text scaling in plots
spy <- 46                  # Number of samples in one cycle (yr)
data(mndvi)                # Load data
t <- as.vector(mndvi$day)  # Days since January 1, 2000
r <- t2rad(t, dpc=dpy)     # Transform days of year to radians
v <- as.vector(mndvi$chl)  # MODIS NDVI over eddy flux tower at Coweeta Hyd. Lab
VX <- vec.x(r,v)           # Avg horizontal vector
VY <- vec.y(r,v)           # Avg vertical vector
vx <- mean(VX, na.rm=TRUE) # Avg horizontal vector
vy <- mean(VY, na.rm=TRUE) # Avg vertical vector
rv <- vec_mag(vx,vy)       # Magnitude (length) of average vector
rv_ang <- vec_ang(vx,vy)   # Angle of the avg vector (phenological median)
avec_ang <- avec_ang(rv_ang)   # Vert opposite of med (avg NDVI min/pheno yr start)
pyStartDay <- rad2d(avec_ang, dpc=dpy) # Cal dy (1-365 beg at Jan 1) marking NDVI min
pyStartIdx <- rad2idx(avec_ang, spc=spy) # Index marking point of least activity
CumNDVI <- sum_cycle(v,pyStartIdx,spy)$cumsum # Cumulative NDVI within each pheno yr
  # Note, output is re-centered, and accumulates starting from the annual NDVI
  # minimum (PyStartIdx) and has c-1 yrs of data due to re-centering.
c.CumNDVI <- c-1           # Number of years in CumNDVI
# Calculate phenologicalc param's for year 2000
cy <- 1                        # The first yr of data (which is 2000 here)
es.idx <- window_idx(CumNDVI,c.CumNDVI,cy,lb,ub)[1] # earlySsn CumNDVI index
ems.idx <- window_idx(CumNDVI,c.CumNDVI,cy,lb,ub)[2] # earlySsn CumNDVI index
ms.idx <- window_idx(CumNDVI,c.CumNDVI,cy,lb,ub)[3] # midSsn CumNDVI index
lms.idx <- window_idx(CumNDVI,c.CumNDVI,cy,lb,ub)[4] # midSsn CumNDVI index
ls.idx <- window_idx(CumNDVI,c.CumNDVI,cy,lb,ub)[5] # lateSsn CumNDVI index
# 3 variables below are days calculated relative to the phenological year, i.e.
# calculated starting from avg NDVI minimum. To get the day relative to the
# calendar year (e.g. 1=Jan 1, 2=Jan 2...) add pyStartDay.
pmets$yr <- cy + 1999
pmets$es <- round(t[es.idx] + pyStartDay)   # Day of calendar yr marking early season milestone
pmets$ems <- round(t[ems.idx] + pyStartDay) # Day of calendar yr marking early-mid season
pmets$ms <- round(t[ms.idx] + pyStartDay)   # Day of calendar yr marking mid season
pmets$lms <- round(t[lms.idx] + pyStartDay) # Day of calendar yr marking late-mid season
pmets$ls <- round(t[ls.idx] + pyStartDay)   # Day of calendar yr marking late season
pmets$s_intv <- pmets$ls - pmets$es # Length of the season from ES to LS
pmets$s_avg <- mean(v[es.idx:ls.idx], na.rm=TRUE) # Mean NDVI for Ssn
pmets$s_sd <- sd(v[es.idx:ls.idx]) # Standard deviation of NDVI for Ssn
pmets$s_mag <- vec_mag(mean(VX[es.idx:ls.idx], na.rm=TRUE),
  mean(VY[es.idx:ls.idx], na.rm=TRUE)) # Magnitude (length) of average vector
    # Magnitude of avg vec between ES & MS thresholds
pmets$ems_mag <- vec_mag(mean(VX[es.idx:(ms.idx-1)], na.rm=TRUE),
  mean(VY[es.idx:(ms.idx-1)], na.rm=TRUE)) # Magnitude of avg vec between ES & MS
pmets$lms_mag <- vec_mag(mean(VX[ms.idx:(ls.idx-1)], na.rm=TRUE),
  mean(VY[ms.idx:(ls.idx-1)], na.rm=TRUE)) # Magnitude of avg vec between MS & LS
print(pmets)

### Example 3, (see Fig. 1 above) Transform into polar coordinates and compare
library(PolarMetrics)
library(plotrix)           # Load plotrix for plotting
dev.new()
dpy <- 365                 # Days/year
c <- 12                    # Num. of years/cycles
data(mndvi)                # Load data
t <- as.vector(mndvi$day)  # Days since January 1, 2000
r <- t2rad(t, dpc=dpy)     # Transform days of year to radians
v <- as.vector(mndvi$wc)   # MODIS NDVI for Willow Creek flux tower, WI
vx <- mean(vec.x(r,v), na.rm=TRUE) # Avg horizontal vector
vy <- mean(vec.y(r,v), na.rm=TRUE) # Avg vertical vector
rv <- vec_mag(vx,vy)       # Magnitude (length) of average vector
rv_ang <- vec_ang(vx,vy)   # Angle of the avg vector (phenological median)
avec_ang <- avec_ang(rv_ang)   # Vert opposite of med (avg NDVI min/pheno yr start)
dev.new(width=12,height=6) # New plot window
par(mfrow=c(1,2))          # Make this a multi-pane plot
cxs <- 1                   # Text scaling in plots
# Time series plot of Willow Creek
plot(2000+(t/dpy), v, pch=20, col='black', xlab='Years',
  ylab='NDVI', main='(a) NDVI Time Series',
  ylim=c(0.35,0.95), cex=cxs, cex.axis=cxs, cex.lab=cxs, cex.main=cxs*1.25)
# Polar plot phenology variables
ndvilabs <- c('','.2','.4','.6','.8','') # labels for radial plot
s.pos <- pi/2 # Radial position to start plotting from
lab.pos <- c(seq(from=0, to=2*pi-(2*pi/12), by=(2*pi)/12))[-4]
rad.labs <- c(month.abb[seq(from=1, to=12)])[-4]
clrs <- colorRampPalette(c('blue3', 'gold'))(length(r)) # Color ramp for plot
# Polar plot of Willow Creek
radial.plot(v,r,clockwise=TRUE,start=s.pos,
  label.pos=lab.pos,labels=rad.labs,radial.labels=ndvilabs,
  rp.type='s',point.symbols=20,point.col=clrs,radial.lim=c(0,1),
	show.radial.grid=FALSE,main='(b) NDVI Polar Plot',
	grid.col='black',grid.unit='NDVI')
	radial.plot(c(0,1),c(0,rv_ang),
		clockwise=TRUE,start=s.pos,rp.type='r',
		lwd=2,line.col='gray45',add=TRUE) # rv_ang, Angle of avg vec
	radial.plot(c(0,1),c(0,avec_ang),
		clockwise=TRUE,start=s.pos,rp.type='s',
		point.symbols='*',cex=cxs*2,
		add=TRUE) # avec_ang, opposite angle of rv_ang
	radial.plot(c(0,rv),c(rv_ang,rv_ang),
		clockwise=TRUE,start=s.pos,rp.type='r',
		lwd=3,line.col='gray45',add=TRUE)		# rv, Magnitude of avg vec
	radial.plot(c(0,rv),c(rv_ang,rv_ang),
		clockwise=TRUE,start=s.pos,rp.type='s',
		point.symbols='*',cex=cxs*2,
		point.col='red',add=TRUE)   # rv, Magnitude of avg vec
	text(sin(avec_ang)*1.15, cos(avec_ang)*1.1,
		'AVec', col='gray45', cex=cxs) # Add text label
	text(sin(rv_ang)*1.2, cos(rv_ang)*1.1,
		'RVec', col='gray45', cex=cxs) # Add text label
	text(0,0.1,'SMag',col='red',cex=cxs) # Add text label
legend('bottomright',c('2000','2011'),col=c('blue3','gold'),
  pch=20,cex=0.7,pt.cex=2,box.lwd=0)
par(xpd=FALSE)
print('Note that 1 year of data is removed, because polar transformation')
print('shifts and centers the data on the phenological year.')


### Example 4, (see Fig. 2 above) Plot cum. NDVI and specify completion milestones
library(PolarMetrics)
dev.new()
dpy <- 365                 # Days/year
c <- 12                    # Num. of years/cycles
spy <- 46                  # Samples/year (or samples/cycle)
data(mndvi)                # Load data
t <- as.vector(mndvi$day)  # Days since January 1, 2000
r <- t2rad(t, dpc=dpy)     # Transform days of year to radians
v <- as.vector(mndvi$wc)   # MODIS NDVI for Willow Creek flux tower, WI
vx <- mean(vec.x(r,v), na.rm=TRUE) # Avg horizontal vector
vy <- mean(vec.y(r,v), na.rm=TRUE) # Avg vertical vector
rv_ang <- vec_ang(vx,vy)   # Angle of the avg vector (phenological median)
avec_ang <- avec_ang(rv_ang)   # Vert opposite of med (avg NDVI min/pheno yr start)
pyStartDay <- rad2d(avec_ang, dpc=dpy) # Cal dy (1-365 beg at Jan 1) marking NDVI min
pyStartIdx <- rad2idx(avec_ang, spc=spy) # Index marking point of least activity
CumNDVI <- sum_cycle(v,pyStartIdx,spy)$cumsum # Cumulative NDVI within each pheno yr
  # Note, output is re-centered, and accumulates starting from the annual NDVI
  # minimum (PyStartIdx) and has c-1 yrs of data due to re-centering.
c.CumNDVI <- c-1           # Number of years in CumNDVI
cy <- 2 # Year within which to examine and plot cumulative NDVI
es <- window_idx(CumNDVI,c.CumNDVI,cy,0.15,0.8)[1] # Indx of CumNDVI marking ES
ms <- window_idx(CumNDVI,c.CumNDVI,cy,0.15,0.8)[3] # Indx of CumNDVI marking MS
ls <- window_idx(CumNDVI,c.CumNDVI,cy,0.15,0.8)[5] # Indx of CumNDVI marking LS
c.CumNDVI <- c-1           # Number of years in CumNDVI
t2 <- seq(pyStartDay,by=8,length.out=47)
y <- c(0,CumNDVI[(spy+1):(2*spy)])
dev.new(width=4,height=4) # New plot window
cxs <- 1                   # Text scaling in plots
plot(t2,y,type='l', xlim=c(15,400),
	xaxt='n', ylim=c(0,max(y)),yaxs='i', yaxt='n', xlab='', ylab='',
	col='gray45', lwd=cxs*2,cex=cxs,cex.axis=cxs,cex.lab=cxs,cex.main=cxs)
mtext('Calculating Timing Metrics', col='black',
	side=3, line=2.3, cex=cxs*1.35) # Main title
y1 <- seq(0,30,by=5)
y2 <- seq(15,85,by=10)
x1 <- seq(0,360,by=30)
x2 <- seq(0,360,by=30)
axis(side=2, labels=y1, at=y1, padj=0.7, col='black', col.axis='black',
	cex.axis=cxs) # Plot left y axis
mtext('Cumulative NDVI', col='black',
	side=2, line=1.9, cex=cxs*0.9)
## x1 axis
axis(side=1, labels=x1, at=x1+pyStartDay, col='black', col.axis='black',
	cex.axis=cxs) # Plot right y axis
## x2 axis
par(tcl = -0.5) # Switch back to outward facing tick marks
axis(side=3, labels=x2, at=x2, col='gray45', col.axis='gray45',
	cex.axis=cxs) # Plot right y axis
par(tcl = 1) # Switch to inward facing tick marks
mtext('Day of Calendar Year',side=3,line=-1,
	col='gray45',cex=cxs*0.9)
mtext('Day of Phenological Year',side=1,line=-1,
	col='black',cex=cxs*0.9)
##
axis(side=4, labels=y2, at=y2/(100/max(y)), col='gray45', col.axis='gray45',
	cex.axis=cxs) # Plot right y axis
mtext('\% Cumulative NDVI',side=4,line=-0.1,
	col='gray45',cex=cxs*0.9)
lines(t2[(es-spy+1):(ls-spy+1)],
	y[(es-spy+1):(ls-spy+1)], col='black', lwd=cxs*3)
#
segments(pyStartDay,0,pyStartDay,max(y),
	col='gray45',lwd=cxs*3, lty=2) # Phen yr offset dahsed line
arrows(pyStartDay,31.5,pyStartDay,max(y),
	col='gray45',lwd=cxs*3) # Phen yr offset arrow head
text(pyStartDay-13,15,srt=90,'Start of Phenological Year',cex=cxs,col='gray45')
#
segments(t2[es-spy+1],y[(es-spy+1)],t2[es-spy+1],30,
	col='gray45',lwd=cxs*3, lty=2) # ES vertical dashed line (calendar yr)
arrows(t2[es-spy+1],30,t2[es-spy+1],32,
	col='gray45',lwd=cxs*3) # ES arrow head (calendar yr)
segments(t2[es-spy+1],y[(es-spy+1)],t2[es-spy+1],4,
	col='black',lwd=cxs*3, lty=2) # ES vertical dashed line (phen yr)
arrows(t2[es-spy+1],4,t2[es-spy+1],2,
	col='black',lwd=cxs*3) # ES arrow head (phen yr)
segments(t2[es-spy+1],y[(es-spy+1)],max(t2),y[(es-spy+1)],
	col='gray45',lwd=cxs*3,lty=2) # ES horizontal dashed line
text(t2[es-spy+1]+23,y[(es-spy+1)]-1,srt=0,'ES',cex=cxs,col='blue3')
#
segments(t2[ms-spy+1],y[(ms-spy+1)],t2[ms-spy+1],30,
	col='gray45',lwd=cxs*3, lty=2) # MS vertical dashed line (calendar yr)
arrows(t2[ms-spy+1],30,t2[ms-spy+1],32,
	col='gray45',lwd=cxs*3) # MS arrow head (calendar yr)
segments(t2[ms-spy+1],y[(ms-spy+1)],t2[ms-spy+1],4,
	col='black',lwd=cxs*3, lty=2) # MS vertical dashed line (phen yr)
arrows(t2[ms-spy+1],4,t2[ms-spy+1],2,
	col='black',lwd=cxs*3) # MS arrow head (phen yr)
segments(t2[ms-spy+1],y[(ms-spy+1)],max(t2),y[(ms-spy+1)],
	col='gray45',lwd=cxs*3,lty=2) # MS horizontal dashed line
text(t2[ms-spy+1]+23,y[(ms-spy+1)]-1,srt=0,'MS',cex=cxs,col='blue3')
#
segments(t2[ls-spy+1],y[(ls-spy+1)],t2[ls-spy+1],30,
	col='gray45',lwd=cxs*3, lty=2) # LS vertical dashed line (calendar yr)
arrows(t2[ls-spy+1],30,t2[ls-spy+1],32,
	col='gray45',lwd=cxs*3) # LS arrow head (calendar yr)
segments(t2[ls-spy+1],y[(ls-spy+1)],t2[ls-spy+1],4,
	col='black',lwd=cxs*3, lty=2) # LS vertical dashed line (phen yr)
arrows(t2[ls-spy+1],4,t2[ls-spy+1],2,
	col='black',lwd=cxs*3) # LS arrow head (phen yr)
segments(t2[ls-spy+1],y[(ls-spy+1)],max(t2),y[(ls-spy+1)],
	col='gray45',lwd=cxs*3,lty=2) # LS horiztontal dashed line
text(t2[ls-spy+1]+23,y[(ls-spy+1)]-1,srt=0,'LS',cex=cxs,col='blue3')
#
par(tcl = -0.5) # Switch back to outward facing tick marks


### Example 5, (see Fig. 3 above) Plot 4 pheno yrs & show phen params of grow.seas.
library(PolarMetrics)
library(plotrix)           # Load plotrix for plotting
dev.new()
lb <- 0.15                 # Lower cumulative sum threshold for window_idx
ub <- 0.8                  # Upper cumulative sum threshold for window_idx
dpy <- 365                 # Days/year
c <- 12                    # Num. of years/cycles
spy <- 46                  # Number of samples in one cycle (yr)
data(mndvi)                # Load data
t <- as.vector(mndvi$day)  # Days since January 1, 2000
r <- t2rad(t, dpc=dpy)     # Transform days of year to radians
v <- as.vector(mndvi$fef)  # MODIS NDVI for pixel covering Fraser Exp. Forest
vx <- mean(vec.x(r,v), na.rm=TRUE) # Avg horizontal vector
vy <- mean(vec.y(r,v), na.rm=TRUE) # Avg vertical vector
rv_ang <- vec_ang(vx,vy)   # Angle of the avg vector (phenological median)
avec_ang <- avec_ang(rv_ang)   # Vert opposite of med (avg NDVI min/pheno yr start)
pyStartIdx <- rad2idx(avec_ang, spc=spy) # Index marking point of least activity
CumNDVI <- sum_cycle(v,pyStartIdx,spy)$cumsum # Cumulative NDVI within each pheno yr
  # Note, output is re-centered, and accumulates starting from the annual NDVI
  # minimum (PyStartIdx) and has c-1 yrs of data due to re-centering.
c.CumNDVI <- c-1           # Number of years in CumNDVI
# Next remove data prior to first NDVI minimum and after last NDVI min
idx=sum_cycle(v,pyStartIdx,spy)$vidx # Indices of v from start of first pheno yr
  # (NDVI min) to end of last complete phenological yr
v1.gs=v[idx]
r1.gs=r[idx]
# Polar plot phenology variables
cxs <- 1                   # Text scaling in plots
ndvilabs <- c('','.2','.4','.6','.8','') # labels for radial plot
s.pos <- pi/2 # Radial position to start plotting from
lab.pos <- c(seq(from=0, to=2*pi-(2*pi/12), by=(2*pi)/12))[-4]
rad.labs <- c(month.abb[seq(from=1, to=12)])[-4]
# Calculate threshold completion milestones for 4 consecutive years
for (I in 1:(c-1)) {
  es.ms.ls=window_idx(
    CumNDVI,c.CumNDVI,I,lb,ub) # Indices for early, mid, late seas.
  idx2=es.ms.ls[1]:es.ms.ls[5] # Indices within this yrs season
  if (I==1) {
    df=cbind(r1.gs[idx2],v1.gs[idx2],I)
  } else {
    df=rbind(df,cbind(r1.gs[idx2],v1.gs[idx2],I))
  }
}
# Combine the calculated polar measures into a data frame
df1=data.frame(df)
colnames(df1)=c('radians','ndvi','yr')
# Plot the 4 years of polar measures on 4 polar plots
dev.new(width=6,height=12)
par(mfrow=c(2,2))
for (I in 1:4) {
  df1.tmp=df1[df1$yr==I,]
  vx2=mean(vec.x(df1.tmp$radians,df1.tmp$ndvi), na.rm=TRUE)
  vy2=mean(vec.y(df1.tmp$radians,df1.tmp$ndvi), na.rm=TRUE)
  rv2=vec_mag(vx2,vy2)
  rv_ang2=vec_ang(vx2,vy2)
  radial.plot(df1.tmp$ndvi,df1.tmp$radians,clockwise=TRUE,start=s.pos,
    label.pos=c(seq(from=0, to=2*pi-(2*pi/12), by=(2*pi)/12)),
    labels=c(month.abb[seq(from=1, to=12)]),
    rp.type='s',point.symbols=20,radial.lim=c(0,1),
    radial.labels=ndvilabs,
    show.radial.grid=FALSE,
    main=paste(I+1999,' Growing Season NDVI',sep=''),
    grid.col='black',point.col='blue3',
    grid.unit='NDVI',cex=cxs*2,cex.axis=cxs,cex.lab=cxs)
  radial.plot(c(0,rv2),c(rv_ang2,rv_ang2),
		clockwise=TRUE,start=s.pos,rp.type='r',
		lwd=3,line.col='red',add=TRUE,cex=cxs) # rv, Magnitude of avg vec
  radial.plot(c(0,rv2),c(rv_ang2,rv_ang2),
		clockwise=TRUE,start=s.pos,rp.type='s',
		point.symbols='*',point.col='red',
		add=TRUE,cex=cxs*2)                        # rv endpoints
  text(-0.15, 0.1,
    'SMag', col='red', cex=cxs) # Vector magnitude
  text(sin(df1.tmp$radians[1])*df1.tmp$ndvi[1]+0.1,
    cos(df1.tmp$radians[1])*df1.tmp$ndvi[1]+0.1,
	  'ES',col='blue3',cex=cxs) # Early season
  text(sin(df1.tmp$radians[nrow(df1.tmp)])*df1.tmp$ndvi[nrow(df1.tmp)]+0.1,
    cos(df1.tmp$radians[nrow(df1.tmp)])*df1.tmp$ndvi[nrow(df1.tmp)]+0.1,
	  'LS',col='blue3',cex=cxs) # Late season
}



}

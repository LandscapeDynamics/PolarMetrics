\name{mndvi}
\alias{mndvi}
\docType{data}
\title{
  MODIS NDVI data from select National Forest sites
}
\description{
  A data set of Normalized Difference Vegetation Index values (NDVI) from MODIS for a select number of locations. NDVI was sampled each 8 days (46 per year) for a total of 552 samples spanning years 2000 through 2011. MODIS sensors have a spatial resolution of approximately 250 square meters, thus each column represents a 250 x 250 meter area.
}
\usage{data("mndvi")}
\format{
  \code{mndvi} is a data frame with 552 rows and 5 columns. Each column corresponds to the variables as follows:
  \describe{
    \item{\code{date}}{Time stamp}
    \item{\code{day}}{a numeric vector of sequential days since January 1, 2000}
    \item{\code{chl}}{a numeric vector of NDVI values for the area surrounding the Coweeta Hydrologic Lab (North Carolina) eddy covariance tower at 83.4252 west longitude and 35.0587 north latitude}
    \item{\code{fef}}{a numeric vector of NDVI values for the area surrounding the Fraser Experimental Forest (Colorado) CO2 tower (run by NCAR) at 105.885417 west longitude, 39.905028 north latitude}
    \item{\code{wc}}{a numeric vector of NDVI values for the area surrounding the Willow Creek (Wisconsin) eddy covariance tower at 90.0798 west longitude, 45.8060 north latitude}
  }
}
%\details{
%}
\source{
  \url{http://dx.doi.org/10.3334/ORNLDAAC/1299}
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(mndvi)
  \dontrun{
  data(mndvi) # Load data
  head(mndvi) # Inspect first few lines of data
  }
}
\keyword{datasets}


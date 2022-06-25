\name{time}

\alias{time}
\alias{time<-}
\alias{time,SpatRaster-method}
\alias{time<-,SpatRaster-method}

  
\title{time of SpatRaster layers}

\description{
Get or set the time of the layers of a SpatRaster. 
}

\usage{
\S4method{time}{SpatRaster}(x, format="")

\S4method{time}{SpatRaster}(x, tstep="")<-value
}

\seealso{\code{\link{depth}}}


\arguments{
  \item{x}{SpatRaster}
  \item{format}{One of "", "seconds" (POSIXlt), "days" (Date), "yearmonths" (expressed as decimal years), or "raw" (either seconds since 1970-01-01-00:00, or year or month number). If "", the returned format is (based on) the format that was used to set the time}
  \item{value}{"Date", "POSIXt", "yearmonth" or numeric}
  \item{tstep}{One of "yearmonths", "years", or "months", ignored when a known time class is used as value}
}

\value{
Date  
}


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   

# Date"
d <- as.Date("2001-05-04") + 0:2
time(s) <- d
time(s)

# POSIX (time stored as seconds)
time(s) <- as.POSIXlt(d)
time(s)

# "raw" time
time(s) <- as.numeric(d)
time(s)
}

\keyword{spatial}


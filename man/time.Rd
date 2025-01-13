\name{time}

\alias{time}
\alias{time,SpatRaster-method}
\alias{time,SpatRasterDataset-method}
\alias{has.time}
\alias{has.time,SpatRaster-method}
\alias{has.time,SpatRasterDataset-method}
\alias{time<-}
\alias{time<-,SpatRaster-method}
\alias{time<-,SpatRasterDataset-method}
\alias{timeInfo}
\alias{timeInfo,SpatRaster-method}
\alias{timeInfo,SpatRasterDataset-method}

  
\title{time of SpatRaster layers}

\description{
Get or set the time of the layers of a SpatRaster. Time can be stored as \code{\link{POSIXlt}} (date and time, with  a resolution of seconds, and a time zone), \code{\link{Date}}, "months", "years", or "yearmonths". 

\code{timeInfo} and \code{has.time} are helper functions to understand what a time data a SpatRaster has.
}

\usage{
\S4method{has.time}{SpatRaster}(x)

\S4method{time}{SpatRaster}(x, format="")

\S4method{time}{SpatRaster}(x, tstep="")<-value

\S4method{timeInfo}{SpatRaster}(x)
}

\seealso{\code{\link{depth}}}


\arguments{
  \item{x}{SpatRaster or SpatRasterDataset}
  \item{format}{One of "", "seconds" (POSIXlt), "days" (Date), "yearmonths" (decimal years), "years", "months". If "", the returned format is (based on) the format that was used to set the time}
  \item{value}{\code{Date}, \code{POSIXt}, \code{yearmon} (defined in package zoo), or numeric}
  \item{tstep}{One of "years", "months", "yearmonths". Used when \code{value} is numeric. Ignored when \code{value} is of type \code{Date}, \code{POSIXt}, or \code{yearmon}}
}

\value{
\code{time}: POSIXlt, Date, or numeric
\code{timeInfo}: data.frame with time step and time zone information (if available)
\code{has.time}: logical
}


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   

# Date"
d <- as.Date("2001-05-04") + 0:2
time(s) <- d
time(s)

# POSIX (date/time with a resolution of seconds)
time(s) <- as.POSIXlt(d)
time(s)

# with time zone
time(s) <- as.POSIXlt(Sys.time(), "America/New_York") + 0:2
time(s)
timeInfo(s)

# years
time(s, tstep="years") <- 2000 + 0:2
s

time(s, tstep="months") <- 1:3
s 
}

\keyword{spatial}


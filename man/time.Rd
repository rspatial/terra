\name{time}

\alias{time}
\alias{time<-}
\alias{time,SpatRaster-method}
\alias{time<-,SpatRaster-method}

  
\title{time of SpatRaster layers}

\description{
Get or set the time of the layers of a SpatRaster. Time can be stored as \code{\link{POSIXlt}} (date and time, with  a resolution of seconds, and a time zone), \code{\link{Date}}, "months", "years", "yearmonths", or "raw". 
}

\usage{
\S4method{time}{SpatRaster}(x, format="")

\S4method{time}{SpatRaster}(x, tstep="")<-value
}

\seealso{\code{\link{depth}}}


\arguments{
  \item{x}{SpatRaster}
  \item{format}{One of "", "seconds" (POSIXlt), "days" (Date), "yearmonths" (decimal years), "years", "months", or "raw". If "", the returned format is (based on) the format that was used to set the time}
  \item{value}{\code{Date}, \code{POSIXt}, \code{yearmonth} (defined in package zoo), or numeric}
  \item{tstep}{One of "years", "months", "yearmonths", or "raw". Used when \code{value} is numeric. Ignored when \code{value} is of type \code{Date}, \code{POSIXt}, or \code{yearmonth}}
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

# POSIX (date/time with a resolution of seconds)
time(s) <- as.POSIXlt(d)
time(s)

# with time zone
time(s) <- as.POSIXlt(Sys.time(), "America/New_York") + 0:2
time(s)

# years
time(s, tstep="years") <- 2000 + 0:2
s

time(s, tstep="months") <- 1:3
s 

# "raw" time
time(s) <- as.numeric(d)
time(s)
}

\keyword{spatial}


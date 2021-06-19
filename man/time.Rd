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
\S4method{time}{SpatRaster}(x)

\S4method{time}{SpatRaster}(x)<-value
}

\seealso{\code{\link{depth}}}


\arguments{
  \item{x}{SpatRaster}
  \item{value}{"Date", "POSIXt", or numeric}
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


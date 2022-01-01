\name{direction}

\alias{direction}
\alias{direction,SpatRaster-method}

\title{Direction}

\description{
The direction (azimuth) to or from the nearest cell that is not \code{NA}. The direction unit is in radians, unless you use argument \code{degrees=TRUE}. 
}

\usage{
\S4method{direction}{SpatRaster}(x, filename='', degrees=FALSE, from=FALSE, ...) 
}

\arguments{
\item{x}{SpatRaster}
\item{filename}{Character. Output filename (optional)}
\item{degrees}{Logical. If \code{FALSE} (the default) the unit of direction is radians.}
\item{from}{Logical. Default is \code{FALSE}. If \code{TRUE}, the direction from (instead of to)  the nearest cell that is not \code{NA} is returned}
\item{...}{Additional arguments as for \code{\link{writeRaster}}}
}



\value{SpatRaster}


\seealso{
\code{\link[terra]{distance}}
}


\examples{
r <- rast(ncol=36,nrow=18)
values(r) <- NA
r[306] <- 1
b <- direction(r) 
plot(b)
}

\keyword{spatial}

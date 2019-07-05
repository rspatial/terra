\name{buffer}

\alias{buffer}
\alias{buffer,SpatRaster-method}


\title{buffer}

\description{
Calculate a buffer around all cells that are not \code{NA} in a SpatRaster, or around the objcts in a SpatVector (t.b.d.)
 
Note that the distance unit of the buffer \code{width} parameter is meters if the SpatRaster CRS is (\code{+proj=longlat}), and in map units (typically also meters) if not. 
}

\usage{
\S4method{buffer}{SpatRaster}(x, width=0, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
\item{x}{SpatRaster}
\item{width}{numeric > 0. Unit is meter if \code{x} has a longitude/latitude CRS, or mapunits in other cases}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{Additional arguments. None implemented}  
}


\value{SpatRaster}

\seealso{
\code{\link{distance}}
}


\examples{
r <- rast(ncol=36,nrow=18)
v <- rep(NA, ncell(r))
v[500] <- 1
values(r) <- v
b <- buffer(r, width=5000000) 
#plot(b)
}

\keyword{spatial}

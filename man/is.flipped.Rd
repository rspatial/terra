\name{is.flipped}

\docType{methods}

\alias{is.flipped}
\alias{is.flipped,SpatRaster-method}
  
  
\title{Is a SpatRaster is flipped}

\description{
Check if a SpatRaster is "flipped" vertically, and may need to be adjusted with \code{\link{flip}} before it can be used. 
}

\usage{
\S4method{is.flipped}{SpatRaster}(x)
}


\arguments{
  \item{x}{SpatRaster}
}  


\value{
logical. One value for each raster data *source*
}


\seealso{
 \code{\link{flip}, \link{is.rotated}}
}

\examples{
r <- rast(nrows=10, ncols=10)
is.flipped(r)
}

\keyword{methods}
\keyword{spatial}


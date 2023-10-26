\name{is.rotated}

\docType{methods}

\alias{is.rotated}
\alias{is.rotated,SpatRaster-method}
  
\title{Check for rotation}

\description{
Check if a SpatRaster is "rotated" and needs to be rectified before it can be used 

See \code{\link{rectify}} 
}

\usage{
\S4method{is.rotated}{SpatRaster}(x)
}


\arguments{
  \item{x}{SpatRaster}
}  


\value{
logical. One value for each raster data *source*
}


\seealso{
 \code{\link{rectify}}
}

\examples{
r <- rast(nrows=10, ncols=10, vals=1:100)
is.rotated(r)
}

\keyword{methods}
\keyword{spatial}


\name{animate}

\docType{methods}

\alias{animate}
\alias{animate,SpatRaster-method}

\title{Animate a SpatRaster}

\description{
Animate (sequentially plot) the layers of a SpatRaster to create a movie.
}

\usage{
\S4method{animate}{SpatRaster}(x, pause=0.25, main, range=NULL, maxcell=50000, n=1, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{pause}{numeric. How long should be the pause be between layers?}  
  \item{main}{title for each layer. If not supplied the z-value is used if available. Otherwise the names are used.}  
  \item{range}{numeric vector of length 2. Range of values to plot, If \code{NULL} the range of all layers is used. If \code{NA} the range of each individual layer is used}
  \item{maxcell}{positive integer. Maximum number of cells to use for the plot. If \code{maxcell < ncell(x)}, \code{spatSample(type="regular")} is used before plotting}
  \item{n}{integer > 0. Number of loops}
  \item{...}{Additional arguments passed to \code{\link{plot}}}
}


\value{
None
}

\seealso{ \code{\link{plot}}}

\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
animate(s, n=1)
}

\keyword{methods}
\keyword{spatial}


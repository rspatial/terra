\name{distValues}

\alias{distValues}
\alias{distValues,SpatRaster,matrix-method}
\alias{distValues,SpatRaster,SpatVector-method}
\alias{distValues,SpatRaster,data.frame-method}


\title{Distance to reference values}

\description{
For each cell, compute the distance in value-space to reference values. References are a matrix or data.frame whose columns match the layer names of \code{x}, or a SpatVector of locations from which values are extracted.
}

\usage{
\S4method{distValues}{SpatRaster,matrix}(x, y, fun="squared", weights=NULL, ...,
		filename="", overwrite=FALSE, wopt=list())

\S4method{distValues}{SpatRaster,data.frame}(x, y, fun="squared", weights=NULL, ...,
		filename="", overwrite=FALSE, wopt=list())

\S4method{distValues}{SpatRaster,SpatVector}(x, y, fun="squared", center=FALSE, scale=FALSE,
		weights=NULL, ..., filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{matrix, data.frame, or SpatVector}
  \item{fun}{character. \code{"abs"} (mean absolute difference) or \code{"squared"} (mean squared difference). Or a function}
  \item{weights}{numeric. Optional weights for the layers of \code{x}}
  \item{center}{logical or numeric. Passed to \code{\link{scale}} (SpatVector method only)}
  \item{scale}{logical or numeric. Passed to \code{\link{scale}} (SpatVector method only)}
  \item{...}{additional arguments passed to \code{fun} (for the built-in functions, \code{na.rm=TRUE})}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{ \code{\link{bestMatch}}, \code{\link{scale}} }

\examples{
r <- rast(system.file("ex/logo.tif", package="terra"))
pts <- vect(cbind(c(25.25, 34.324), c(54.577, 46.489)))
x <- scale(r)
d <- distValues(x, pts)
}

\keyword{spatial}

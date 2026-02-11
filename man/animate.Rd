\name{animate}

\docType{methods}

\alias{animate}
\alias{animate,SpatRaster-method}
\alias{animate,SpatVector-method}

\title{Animate a map}

\description{
Animate (sequentially plot) the layers of a SpatRaster, or the variables or geometries of a SpatVector, to create a movie.
}

\usage{
\S4method{animate}{SpatRaster}(x, pause=0.25, main, range=NULL, maxcell=50000, n=1, ...)

\S4method{animate}{SpatVector}(x, pause=0.25, main="", n=1, vars=NULL, range=NULL, add=NULL, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{pause}{numeric. How long should the pause be between layers?}  
  \item{main}{title for each layer. For SpatRaster, if not supplied, the z-value is used if available. Otherwise the names are used.}  
  \item{range}{numeric vector of length 2. Range of values to plot, If \code{NULL} the range of all layers is used for rasters, or all variables for vectors if they are all numeric. If \code{NA} the range of each individual layer is used}
  \item{maxcell}{positive integer. Maximum number of cells to use for the plot. If \code{maxcell < ncell(x)}, \code{spatSample(type="regular")} is used before plotting}
  \item{n}{integer > 0. Number of plotting loops}
  \item{vars}{numeric or character to indicate the variables to animate. If this is NULL, the geometries are animated instead}
  \item{add}{logical. Add the geometries to the current map? When looping over variables: \code{NULL} is equivalent to \code{TRUE}. When looping over geometries: if \code{TRUE}, add all geometries to the current plot. If \code{NULL}, \code{add} is set to \code{FALSE} for the first geometry and \code{TRUE} for the remaining ones.}
  \item{...}{additional arguments passed to \code{\link{plot}}}
}


\value{
None
}

\author{Márcia Barbosa, Robert J. Hijmans}

\seealso{ \code{\link{plot}}}

\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
animate(s, n=1)

v <- vect(system.file("ex/lux.shp", package="terra"))
animate(v, n=2)
animate(v, n=1, vars=names(v))

# you can save an animation to file like this
# animation::saveGIF(terra::animate(v, n=1), "animation.gif")
}

\keyword{methods}
\keyword{spatial}


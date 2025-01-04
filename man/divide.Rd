\name{divide}

\alias{divide}
\alias{divide,SpatRaster-method}
\alias{divide,SpatVector-method}

\title{
Subdivide a raster or polygons
}

\description{
Divide a \code{SpatRaster} into \code{n} parts with approximately the same sum of weights (cell values). 

Divides a \code{SpatVector} of polygons into \code{n} compact and approximately equal area parts. The results are not deterministic so you should use set.seed to be able to reproduce your results. If you get a warning about non-convergence, you can increase the number of iterations used with additional argument iter.max
}

\usage{
\S4method{divide}{SpatRaster}(x, n=2, start="ns", as.raster=FALSE, na.rm=TRUE)

\S4method{divide}{SpatVector}(x, n=5, w=NULL, alpha=1, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatVector of polygons}
  \item{n}{numeric. Can be a single positive integer to indicate the number of parts (SpatVector) or the number of splits (SpatRaster). 
  
  If \code{x} is a SpatRaster, it can also be a vector with values -2, -1, 1, or 2. Where 1 means one split and 2 means two splits, and the negative sign indicates an East-West (vertical) split as opposed to a North-South split. 
  
  If \code{x} is a SpatVector it can be a list with at least one of these elements: \code{horizontal} and \code{vertical} that specify the proportions of the area that splits should cover. This can either be a single fraction such as 1/3, or a sequence of fractions in ascending order such as \code{c(1/4, 1/2, 1)}}
  
  \item{start}{character. To indicate the initial direction of splitting the raster. "ns" for North-South (horizontal) or "ew" for East-West (vertical)}
  \item{as.raster}{logical. If \code{FALSE} a SpatVector is returned. If \code{FALSE}, a SpatRaster is returned. If \code{NA} a list with a SpatRaster and a SpatVector is returned}
  \item{na.rm}{logical. If \code{TRUE} cells in \code{x} that are \code{NA} are not included in the output}  

  \item{w}{SpatRaster with, for example, environmental data}
  \item{alpha}{numeric. One or two numbers that act as weights for the x and y coordinates}
  \item{...}{additional arguments such as \code{iter.max} passed on to \code{\link{kmeans}}}

}


\value{
SpatVector or SpatRaster, or a list with both
}


\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
x <- divide(r, 3)
# plot(r); lines(x)


f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
d <- divide(v, 3)
dv <- divide(v, list(h=.5))
}

\keyword{spatial}


\name{divide}

\alias{divide}
\alias{divide,SpatRaster-method}
\alias{divide,SpatVector-method}

\title{
Subdivide a raster or polygons
}

\description{
Divide a \code{SpatRaster} by splitting it \code{n} times. By default the resulting parts have approximately the same sum of cell values (\code{split="values"}). Use \code{split="notNA"} for the same number of non-\code{NA} cells, or \code{split="cells"} to ignore values and \code{NA}s and split the grid by cell count.

Divide a \code{SpatVector} of polygons into \code{n} compact and approximately equal-area parts. The results are not deterministic, so you should use \code{\link{set.seed}} to be able to reproduce your results. If you get a warning about non-convergence, you can increase the number of iterations used with additional argument \code{iter.max}
}

\usage{
\S4method{divide}{SpatRaster}(x, n=2, start="ns", as.raster=FALSE, na.rm=TRUE, split="values")

\S4method{divide}{SpatVector}(x, n=5, w=NULL, alpha=1, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatVector of polygons}
  \item{n}{numeric. 
  
  For a SpatVector, a single positive integer: the number of parts.
  
  For a SpatRaster, a single positive integer is the number of split rounds. Each round cuts every current piece in two, alternating direction (see \code{start}). So \code{n=1} gives 2 parts, \code{n=2} gives 4, \code{n=3} gives 8.
  
  For a SpatRaster, \code{n} can instead be a vector that lists the cuts, in order. Each value is applied to every piece from the previous step. The sign is the direction: positive = North-South (a horizontal cut), negative = East-West (a vertical cut). The absolute value is the number of cuts in that step, so the piece is divided into \code{abs(n)+1} parts (\code{1} = halves, \code{2} = thirds, \code{3} = quarters, ...). For example \code{c(1, -1)} is a North-South cut into 2, then an East-West cut of each half (4 parts); \code{c(2, -1, 1)} is a North-South cut into 3, then East-West into 2, then North-South into 2 (12 parts).
  
  For a SpatVector, \code{n} can be a list with at least one of these elements: \code{horizontal} and \code{vertical} that specify the proportions of the area that splits should cover. This can either be a single fraction such as 1/3, or a sequence of fractions in ascending order such as \code{c(1/4, 1/2, 1)}}
  
  \item{start}{character. Initial split direction when \code{n} is a single number. \code{"ns"} for North-South (horizontal cuts) or \code{"ew"} for East-West (vertical cuts). Ignored if \code{n} is a vector}
  \item{as.raster}{logical. If \code{FALSE} a SpatVector is returned. If \code{TRUE}, a SpatRaster is returned. If \code{NA} a list with a SpatRaster and a SpatVector is returned}
  \item{na.rm}{logical. If \code{TRUE} cells in \code{x} that are \code{NA} are not included in the output}
  \item{split}{character. For SpatRaster only. \code{"values"} (the default) balances the sum of cell values in each part. \code{"notNA"} balances the number of non-\code{NA} cells. \code{"cells"} ignores values and \code{NA}s and splits by the number of cells}  

  \item{w}{SpatRaster with, for example, environmental data}
  \item{alpha}{numeric. One or two numbers that act as weights for the x and y coordinates}
  \item{...}{additional arguments such as \code{iter.max} passed on to \code{\link{kmeans}}}
}

\seealso{\code{\link{thresh}}; \code{\link{makeTiles}}; \code{\link{split}}; \code{\link{segregate}}}

\value{SpatVector or SpatRaster, or a list with both}


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


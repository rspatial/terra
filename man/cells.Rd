\name{cells}

\docType{methods}

\alias{cells}
\alias{cells,SpatRaster,missing-method}
\alias{cells,SpatRaster,numeric-method}
\alias{cells,SpatRaster,SpatVector-method}
\alias{cells,SpatRaster,SpatExtent-method}


\title{Get cell numbers}

\description{ 
Get the cell numbers covered by a SpatVector or SpatExtent. Or that match values in a vector; or all non \code{NA} values.
}

\usage{
\S4method{cells}{SpatRaster,missing}(x, y)

\S4method{cells}{SpatRaster,numeric}(x, y, pairs=FALSE)

\S4method{cells}{SpatRaster,SpatVector}(x, y, method="simple", weights=FALSE, exact=FALSE, 
				touches=is.lines(y), small=TRUE)

\S4method{cells}{SpatRaster,SpatExtent}(x, y)
}


\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatVector, SpatExtent, 2-column matrix representing points, numeric representing values to match, or missing}
  \item{method}{character. Method for getting cell numbers for points. The default is "simple", the alternative is "bilinear". If it is "bilinear", the four nearest cells and their weights are returned}
  \item{weights}{logical. If \code{TRUE} and \code{y} has polygons, the approximate fraction of each cell that is covered is returned as well}
  \item{pairs}{logical. If \code{TRUE} the cell values matched area also returned}
  \item{exact}{logical. If \code{TRUE} and \code{y} has polygons, the exact fraction of each cell that is covered is returned as well}
  \item{touches}{logical. If \code{TRUE}, values for all cells touched by lines or polygons are extracted, not just those on the line render path, or whose center point is within the polygon. Not relevant for points}
  \item{small}{logical. If \code{TRUE}, values for all cells in touched polygons are extracted if none of the cells center points is within the polygon; even if \code{touches=FALSE}}
}

\value{
numeric vector or matrix
}


\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
r[c(1:25, 31:100)] <- NA
r <- ifel(r > 28, r + 10, r)

# all cell numbers of cells that are not NA
cells(r)

# cell numbers that match values
x <- cells(r, c(28,38))
x$lyr.1

# cells for points
m <- cbind(x=c(0,10,-30), y=c(40,-10,20))
cellFromXY(r, m)

v <- vect(m)
cells(r, v)
cells(r, v, method="bilinear")

# cells for polygons
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
r <- rast(v)
cv <- cells(r, v) 
}


\keyword{methods}
\keyword{spatial}


\name{as.spatvector}

\docType{methods}

\alias{as.polygons}
\alias{as.lines}
\alias{as.points}

\alias{as.polygons,SpatRaster-method}
\alias{as.lines,SpatRaster-method}
\alias{as.points,SpatRaster-method}
\alias{as.polygons,SpatVector-method}
\alias{as.lines,SpatVector-method}
\alias{as.points,SpatVector-method}
\alias{as.polygons,SpatExtent-method}
\alias{as.lines,SpatExtent-method}
\alias{as.points,SpatExtent-method}

\title{Conversion to a SpatVector, or to another SpatVector type}

\description{
Conversion of a SpatRaster or SpatExtent to a SpatVector of points, lines, or polygons; 

And conversion of a SpatVector to another SpatVector type.
}

\usage{
\S4method{as.polygons}{SpatRaster}(x, trunc=TRUE, dissolve=TRUE, values=TRUE,
		na.rm=TRUE, na.all=FALSE, extent=FALSE)

\S4method{as.lines}{SpatRaster}(x)

\S4method{as.points}{SpatRaster}(x, values=TRUE, na.rm=TRUE, na.all=FALSE)

\S4method{as.polygons}{SpatVector}(x, extent=FALSE)

\S4method{as.lines}{SpatVector}(x)

\S4method{as.points}{SpatVector}(x, multi=FALSE, skiplast=TRUE)

\S4method{as.polygons}{SpatExtent}(x, crs="")

\S4method{as.lines}{SpatExtent}(x, crs="")

\S4method{as.points}{SpatExtent}(x, crs="")
}


\arguments{
\item{x}{SpatRaster or SpatVector}
\item{trunc}{logical; truncate values to integers. Cells with the same value are merged. Therefore, if \code{trunc=FALSE} the object returned can have many cells and can be very large}
\item{dissolve}{logical; combine cells with the same values? If \code{TRUE} only the first layer in \code{x} is processed}
\item{values}{logical; include cell values as attributes?}
\item{multi}{logical. If \code{TRUE} a multi-point geometry is returned}
\item{skiplast}{logical. If \code{TRUE} the last point of a polygon (which is the same as the first point) is not included}
\item{extent}{logical. if \code{TRUE}, a polygon for the extent of the SpatRaster or SpatVector is returned. If \code{x} is a SpatRaster, the polygon has vertices for each row and column, not just the four corners of the raster. This can be useful for more precise projection. If that is not required, it is more efficient to get the extent represented by only the four corners with \code{as.polygons(ext(x), crs=crs(x))}}
\item{na.rm}{logical. If \code{TRUE} cells that are \code{NA} are ignored}
\item{na.all}{logical. If \code{TRUE} cells are only ignored if \code{na.rm=TRUE} and their value is \code{NA} for \bold{all} layers instead of for \code{any} layer}
\item{crs}{character. The coordinate reference system (see \code{\link{crs}})}
}

\value{
SpatVector
}


\examples{
r <- rast(ncols=2, nrows=2)
values(r) <- 1:ncell(r)

as.points(r)
as.lines(ext(r), crs=crs(r))

if  (gdal() >= "3.0.0") {
	p <- as.polygons(r)
	p
	as.lines(p)
	as.points(p)
}
}

\keyword{spatial}


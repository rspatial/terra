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

And conversion of a SpatVector to a another SpatVector type.
}

\usage{
\S4method{as.polygons}{SpatRaster}(x, trunc=TRUE, dissolve=TRUE, values=TRUE, na.rm=TRUE, extent=FALSE)

\S4method{as.lines}{SpatRaster}(x)

\S4method{as.points}{SpatRaster}(x, values=TRUE, na.rm=TRUE)

\S4method{as.polygons}{SpatVector}(x)

\S4method{as.lines}{SpatVector}(x)

\S4method{as.points}{SpatVector}(x, multi=FALSE)

\S4method{as.polygons}{SpatExtent}(x, crs="")

\S4method{as.lines}{SpatExtent}(x, crs="")

\S4method{as.points}{SpatExtent}(x, crs="")
}


\arguments{
\item{x}{SpatRaster or SpatVector}
\item{trunc}{logical; truncate values to integers. Cels with the same value are merged. Therefore, if \code{trunc=FALSE} the object returned can be very large}
\item{dissolve}{logical; combine cells with the same values? If \code{TRUE} only the first layer in \code{x} is processed}
\item{values}{logical; include cell values as attributes? If \code{FALSE} the cells are not dissolved and the object returned can be very large}
\item{multi}{logical. If \code{TRUE} a multipoint geometry is returned}
\item{extent}{logical. if \code{TRUE}, a polygon for the extent of the SpatRaster is returned. It has vertices for each grid cell, not just the four corners of the raster. This can be useful for more precise projection. In other cases it is better to do \code{as.polygons(ext(x))} to get a much smaller object returned that covers the same extent}
\item{na.rm}{logical. If \code{TRUE} cells that are \code{NA} are ignored}
\item{crs}{character. The coordinate reference system (see \code{\link{crs}}}
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


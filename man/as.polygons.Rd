\name{as.polygons}

\docType{methods}

\alias{as.polygons}

\alias{as.polygons,SpatRaster-method}
\alias{as.polygons,SpatVector-method}
\alias{as.polygons,SpatExtent-method}

\title{Conversion to a SpatVector of polygons}

\description{
Conversion of a SpatRaster, SpatVector or SpatExtent to a SpatVector of polygons. 
}

\usage{
\S4method{as.polygons}{SpatRaster}(x, round=TRUE, aggregate=TRUE, values=TRUE,
		na.rm=TRUE, na.all=FALSE, extent=FALSE, digits=0, ...)

\S4method{as.polygons}{SpatVector}(x, extent=FALSE)

\S4method{as.polygons}{SpatExtent}(x, crs="")
}


\arguments{
\item{x}{SpatRaster, SpatVector or SpatExtent}
\item{round}{logical; If \code{TRUE} and \code{aggregate=TRUE}, values are rounded before aggregation. If this value is \code{FALSE} the SpatVector returned can have very many polygons and can be very large}
\item{aggregate}{logical; combine cells with the same values? If \code{TRUE} only the first layer in \code{x} is processed}
\item{values}{logical; include cell values as attributes?}
\item{extent}{logical. if \code{TRUE}, a polygon for the extent of the SpatRaster or SpatVector is returned. If \code{x} is a SpatRaster, the polygon has vertices for each row and column, not just the four corners of the raster. This can be useful for more precise projection. If that is not required, it is more efficient to get the extent represented by only the four corners with \code{as.polygons(ext(x), crs=crs(x))}}
\item{na.rm}{logical. If \code{TRUE} cells that are \code{NA} are ignored}
\item{na.all}{logical. If \code{TRUE} cells are only ignored if \code{na.rm=TRUE} and their value is \code{NA} for \bold{all} layers instead of for \code{any} layer}
\item{digits}{integer. The number of digits for rounding (if \code{round=TRUE})}
\item{crs}{character. The coordinate reference system (see \code{\link{crs}})}
\item{...}{additional arguments. For backward compatibility. Will be removed in the future}
}

\value{
SpatVector
}

\seealso{
\code{\link{as.lines}}, \code{\link{as.points}}
}

\examples{
r <- rast(ncols=2, nrows=2)
values(r) <- 1:ncell(r)

if  (gdal() >= "3.0.0") {
	p <- as.polygons(r)
	p
}
}

\keyword{spatial}


\name{as.points}

\docType{methods}

\alias{as.points}

\alias{as.points,SpatRaster-method}
\alias{as.points,SpatVector-method}
\alias{as.points,SpatExtent-method}

\title{Conversion to a SpatVector of points}

\description{
Conversion of a SpatRaster, SpatVector or SpatExtent to a SpatVector of points.
}

\usage{
\S4method{as.points}{SpatRaster}(x, values=TRUE, na.rm=TRUE, na.all=FALSE)

\S4method{as.points}{SpatVector}(x, multi=FALSE, skiplast=TRUE)

\S4method{as.points}{SpatExtent}(x, crs="")
}


\arguments{
\item{x}{SpatRaster, SpatVector or SpatExtent}
\item{values}{logical; include cell values as attributes?}
\item{multi}{logical. If \code{TRUE} a multi-point geometry is returned}
\item{skiplast}{logical. If \code{TRUE} the last point of a polygon (which is the same as the first point) is not included}
\item{na.rm}{logical. If \code{TRUE} cells that are \code{NA} are ignored}
\item{na.all}{logical. If \code{TRUE} cells are only ignored if \code{na.rm=TRUE} and their value is \code{NA} for \bold{all} layers instead of for \code{any} layer}
\item{crs}{character. The coordinate reference system (see \code{\link{crs}})}
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

as.points(r)

if  (gdal() >= "3.0.0") {
	p <- as.polygons(r)
	as.points(p)
}
}

\keyword{spatial}


\name{as.lines}

\docType{methods}

\alias{as.lines}

\alias{as.lines,SpatRaster-method}
\alias{as.lines,SpatVector-method}
\alias{as.lines,SpatExtent-method}
\alias{as.lines,matrix-method}

\title{Conversion to a SpatVector of lines}

\description{
Conversion of a SpatRaster, SpatVector or SpatExtent to a SpatVector of lines. 
}

\usage{
\S4method{as.lines}{SpatRaster}(x)

\S4method{as.lines}{SpatVector}(x)

\S4method{as.lines}{SpatExtent}(x, crs="")

\S4method{as.lines}{matrix}(x, crs="")
}


\arguments{
\item{x}{SpatRaster, SpatVector, SpatExtent or matrix. If \code{x} is a matrix it should have two columns for a single line, or four columns, where each row has the start and end coordinates (x, y) for lines}
\item{crs}{character. The coordinate reference system (see \code{\link{crs}})}
}

\value{
SpatVector
}

\seealso{
\code{\link{as.points}}, \code{\link{as.polygons}}
}


\examples{
r <- rast(ncols=2, nrows=2)
values(r) <- 1:ncell(r)

as.lines(r)

as.lines(ext(r), crs=crs(r))

if  (gdal() >= "3.0.0") {
	p <- as.polygons(r)
	as.lines(p)
}


## with a matrix
s <- cbind(1:5, 1:5)
e <- cbind(1:5, 0)

as.lines(s)
as.lines(cbind(s, e), "+proj=longlat")

}

\keyword{spatial}


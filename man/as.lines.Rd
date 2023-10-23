\name{as.lines}

\docType{methods}

\alias{as.lines}

\alias{as.lines,SpatRaster-method}
\alias{as.lines,SpatVector-method}
\alias{as.lines,SpatExtent-method}

\title{Conversion to a SpatVector of lines}

\description{
Conversion of a SpatRaster, SpatVector or SpatExtent to a SpatVector of lines. 
}

\usage{
\S4method{as.lines}{SpatRaster}(x)

\S4method{as.lines}{SpatVector}(x)

\S4method{as.lines}{SpatExtent}(x, crs="")
}


\arguments{
\item{x}{SpatRaster, SpatVector or SpatExtent}
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
}

\keyword{spatial}


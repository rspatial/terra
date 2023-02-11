\name{geomtype}

\docType{methods}

\alias{geomtype}
\alias{geomtype,SpatVector-method}
\alias{geomtype,SpatVectorProxy-method}
\alias{geomtype,Spatial-method}
\alias{datatype}
\alias{datatype,SpatVector-method}
\alias{datatype,SpatRaster-method}
\alias{is.points}
\alias{is.lines}
\alias{is.polygons}
\alias{is.points,SpatVector-method}
\alias{is.lines,SpatVector-method}
\alias{is.polygons,SpatVector-method}

\title{Geometry type of a SpatVector}

\description{
Get the geometry type (points, lines, or polygons) of a SpatVector or the data types of the fields (attributes, variables) of a SpatVector.
} 

\usage{
\S4method{geomtype}{SpatVector}(x)

\S4method{datatype}{SpatVector}(x)

\S4method{is.points}{SpatVector}(x)

\S4method{is.lines}{SpatVector}(x)

\S4method{is.polygons}{SpatVector}(x)

\S4method{datatype}{SpatRaster}(x, bylyr=TRUE)
}

\arguments{
\item{x}{SpatVector or SpatRaster}

\item{bylyr}{logical. If \code{TRUE} a value is returned for each layer. Otherwise, a value is returned for each data source (such as a file)}
}

\value{character}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)

geomtype(v)
is.polygons(v)
is.lines(v)
is.points(v)

names(v)
datatype(v)
}

\keyword{methods}
\keyword{spatial}


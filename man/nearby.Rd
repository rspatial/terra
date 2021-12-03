\name{nearest}

\docType{methods}

\alias{nearby}
\alias{nearby,SpatVector-method}
\alias{nearest}
\alias{nearest,SpatVector-method}

\title{nearby geometries}

\description{
Identify geometries that are near to each other. Either get the index of all geometries within a certain distance, or the k nearest neighbors, or (with \code{nearest}) get the nearest points between two geometries.
}

\usage{
\S4method{nearby}{SpatVector}(x, y=NULL, distance=0, k=1, centroids=TRUE, symmetrical=TRUE)

\S4method{nearest}{SpatVector}(x, y, pairs=FALSE, centroids=TRUE, lines=FALSE)
}

\arguments{
  \item{x}{SpatVector}
  \item{y}{SpatVector or NULL}
  \item{distance}{numeric. maximum distance}
  \item{k}{positive integer. number of neighbors. Ignored if \code{distance > 0}}
  \item{centroids}{logical. Should the centroids of polygons be used?}
  \item{symmetrical}{logical. If \code{TRUE}, a near pair is only included once. That is, if geometry 1 is near to geometry 3, the implied nearness between 3 and 1 is not reported. Ignored if \code{k} neighbors are returned}
  \item{pairs}{logical. If \code{TRUE} pairwise nearest points are returned (only relevant when using at least one SpatVector of lines or polygons}
  \item{lines}{logical. If \code{TRUE} lines between the nearest points instead of (the nearest) points }
}

\seealso{\code{\link{relate}}, \code{\link{adjacent}}}

\value{
matrix
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
nearby(v, distance=12000)
}

\keyword{spatial}

\name{near}

\docType{methods}

\alias{near}
\alias{near,SpatVector-method}

\title{nearby geometries}

\description{
Identify geometries that are near to each other. Either get the index of all geometries within a cetain distance, or get the k nearest neighbors.
}

\usage{
\S4method{near}{SpatVector}(x, distance=0, k=1, centroids=TRUE, symmetrical=TRUE, ...)
}

\arguments{
  \item{x}{SpatVector}
  \item{distance}{numeric. maximum distance}
  \item{k}{postive integer. number of neighbors. Ignored of distance > 0}
  \item{centroids}{logical. Should the centroids of polygons be used?}
  \item{symmetrical}{logical. If \code{TRUE}, a near pair is only included once. That is, if geometry 1 is near to geometry 3, the implied nearness between 3 and 1 is not reported}
  \item{...}{additional arguments. None implemented}
}

\seealso{\code{\link{relate}}, \code{\link{adjacent}}}

\value{
matrix
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
near(v, distance=12000)
}

\keyword{spatial}

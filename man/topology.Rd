\name{topology}

\docType{methods}

\alias{topology}

\alias{makeNodes}
\alias{mergeLines}
\alias{removeDupNodes}
\alias{emptyGeoms}
\alias{snap}

\alias{makeNodes,SpatVector-method}
\alias{mergeLines,SpatVector-method}
\alias{snap,SpatVector-method}
\alias{removeDupNodes,SpatVector-method}
\alias{emptyGeoms,SpatVector-method}

\title{ 
Vector topology methods
}

\description{

\code{makeNodes} create nodes on lines

\code{mergeLines} connect lines to form polygons

\code{removeDupNodes} removes duplicate nodes in geometries and optionally rounds the coordinates

\code{emptyGeoms} returns the indices of empty (null) geometries. \code{\link{is.na}} also checks if any of the coordinates is \code{NA}.

\code{snap} makes boundaries of geometries identical if they are very close to each other.
}

\usage{
\S4method{mergeLines}{SpatVector}(x)
\S4method{snap}{SpatVector}(x, y=NULL, tolerance)
\S4method{removeDupNodes}{SpatVector}(x, digits = -1)
\S4method{makeNodes}{SpatVector}(x)
}

\arguments{
  \item{x}{SpatVector of lines or polygons}
  \item{y}{SpatVector of lines or polygons to snap to. If \code{NULL} snapping is to the other geometries in \code{x}}
  \item{tolerance}{numeric. Snapping tolerance (distance between geometries)}
  \item{digits}{numeric. Number of digits used in rounding. Ignored if < 0}
} 

\value{
SpatVector
}


\seealso{\code{\link{sharedPaths}}, \code{\link{gaps}}, \code{\link{simplifyGeom}}, \code{\link{forceCCW}}}




\examples{
p1 <- as.polygons(ext(0,1,0,1))
p2 <- as.polygons(ext(1.1,2,0,1))

p <- rbind(p1, p2)

y <- snap(p, tol=.15)
plot(p, lwd=3, col="light gray")
lines(y, col="red", lwd=2)
}

\keyword{methods}
\keyword{spatial}


\name{thin}

\docType{methods}

\alias{thin}
\alias{thin,SpatVector-method}

\alias{thinNodes}
\alias{thinNodes,SpatVector-method}

\title{Subset geometries by minimum distance}

\description{
\code{thin} return a subset of geometries such that all remaining geometries are at least \code{d} apart from each other. This uses greedy spatial thinning: the first geometry is always kept, and subsequent geometries are only kept if they are at least \code{d} away from all previously kept geometries.

\code{thinNodes} reduces the number of nodes to represent line or polygon geometries
}

\usage{
\S4method{thin}{SpatVector}(x, d, unit="m")

\S4method{thinNodes}{SpatVector}(x, d, unit="m", makeValid=TRUE)
}

\arguments{
  \item{x}{SpatVector}
  \item{d}{positive numeric. The minimum distance between geometries or nodes}
  \item{unit}{character. \code{"m"} (meter, the default) or \code{"km"} (kilometer)}
  \item{makeValid}{logical. If \code{TRUE} an attempt is made to assure that thinned polygon geometries are valid (e.g., not self intersecting)}
}

\value{
SpatVector
}

\seealso{\code{\link{densify}}, \code{\link{distance}}, \code{\link{nearby}}, \code{\link{nearest}}}

\examples{
p <- vect(rbind(c(0,0), c(0.5,0.5), c(0.6, 0.6), c(5,5), c(5.5, 5.5), c(50, 50)), crs="+proj=longlat +datum=WGS84")
values(p) <- data.frame(id=1:6)
nrow(p)

# points that are at least 200 km apart
p2 <- thin(p, 200000)
nrow(p2)

# same, using km
p3 <- thin(p, 200, unit="km")
nrow(p3)
}

\keyword{spatial}

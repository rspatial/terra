\name{intersect}

\docType{methods}

\alias{intersect}

\alias{intersect,SpatVector,SpatVector-method}
\alias{intersect,SpatVector,SpatExtent-method}
\alias{intersect,SpatExtent,SpatVector-method}
\alias{intersect,SpatExtent,SpatExtent-method}

\title{ 
Intersection
}

\description{
Intersect the geometries of two SpatVectors.

Intersecting points with points uses the extent of \code{y} to get the intersection. Intersecting of points and lines is not supported because of numerical inaccuracies with that. You can use \code{\link{buffer}}, to create polygons from lines and use these with intersect. 

See \code{\link{crop}} for intersection of a SpatRaster.
}

\usage{
\S4method{intersect}{SpatVector,SpatVector}(x, y)

\S4method{intersect}{SpatVector,SpatExtent}(x, y)

\S4method{intersect}{SpatExtent,SpatVector}(x, y)

\S4method{intersect}{SpatExtent,SpatExtent}(x, y)
}

\arguments{
  \item{x}{SpatVector or SpatExtent}
  \item{y}{SpatVector or SpatExtent}
} 

\value{
Same as \code{x}
}

\seealso{
\code{\link{union}}, \code{\link{crop}}, \code{\link{relate}}
}

\examples{
e1 <- ext(-10, 10, -20, 20)
e2 <- ext(0, 20, -40, 5)
intersect(e1, e2)

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
e <- ext(5.6, 6, 49.55, 49.7)
x <- intersect(v, e)

p <- vect(c("POLYGON ((5.8 49.8, 6 49.9, 6.15 49.8, 6 49.6, 5.8 49.8))", 
"POLYGON ((6.3 49.9, 6.2 49.7, 6.3 49.6, 6.5 49.8, 6.3 49.9))"), crs=crs(v))
values(p) <- data.frame(pid=1:2, area=expanse(p))

y <- intersect(v, p)
}


\keyword{methods}
\keyword{spatial}


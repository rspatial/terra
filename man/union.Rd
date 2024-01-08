\name{union}

\docType{methods}

\alias{union}
\alias{union,SpatExtent,SpatExtent-method}
\alias{union,SpatVector,SpatVector-method}
\alias{union,SpatVector,SpatExtent-method}
\alias{union,SpatVector,missing-method}

\title{ 
Union SpatVector or SpatExtent objects
}

\description{

Overlapping polygons (between, not within, objects) are intersected. Union for lines and points simply combines the two data sets; without any geometric intersections. This is equivalent to \code{\link{c}}. Attributes are joined. See \code{\link{c}} if you want to combine polygons without intersection. 

If \code{x} and \code{y} have a different geometry type, a SpatVectorCollection is returned.

If a single SpatVector is supplied, overlapping polygons are intersected. Original attributes are lost. New attributes allow for determining how many, and which, polygons overlapped. 

SpatExtent: Objects are combined into their union; this is equivalent to \code{+}.  
}

\usage{
\S4method{union}{SpatVector,SpatVector}(x, y)

\S4method{union}{SpatVector,missing}(x, y)

\S4method{union}{SpatExtent,SpatExtent}(x, y)
}

\arguments{
  \item{x}{SpatVector or SpatExtent}
  \item{y}{Same as \code{x} or missing}
} 

\value{
SpatVector or SpatExtent
}

\seealso{
\code{\link[terra]{intersect}}

\code{\link{merge}} and \code{\link{mosaic}} to union SpatRasters.

\code{\link{crop}} and \code{\link{extend}} for the union of SpatRaster and SpatExtent.

\code{\link[terra]{merge}} for merging a data.frame with attributes of a SpatVector.
  
\code{\link[terra]{aggregate}} to dissolve SpatVector objects. 
}


\examples{
e1 <- ext(-10, 10, -20, 20)
e2 <- ext(0, 20, -40, 5)
union(e1, e2)

#SpatVector
v <- vect(system.file("ex/lux.shp", package="terra"))
v <- v[,3:4]
p <- vect(c("POLYGON ((5.8 49.8, 6 49.9, 6.15 49.8, 6 49.65, 5.8 49.8))", 
"POLYGON ((6.3 49.9, 6.2 49.7, 6.3 49.6, 6.5 49.8, 6.3 49.9))"), crs=crs(v))
values(p) <- data.frame(pid=1:2, value=expanse(p))
u <- union(v, p)
plot(u, "pid")

b <- buffer(v, 1000)

u <- union(b)
u$sum <- rowSums(as.data.frame(u))
plot(u, "sum")
}

\keyword{methods}
\keyword{spatial}


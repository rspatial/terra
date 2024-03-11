\name{relate}

\docType{methods}

\alias{relate}
\alias{is.related}

\alias{relate,SpatVector,SpatVector-method}
\alias{relate,SpatVector,missing-method}
\alias{relate,SpatVector,SpatExtent-method}
\alias{relate,SpatExtent,SpatVector-method}
\alias{relate,SpatExtent,SpatExtent-method}
\alias{relate,SpatVector,SpatRaster-method}
\alias{relate,SpatRaster,SpatVector-method}
\alias{relate,SpatExtent,SpatRaster-method}
\alias{relate,SpatRaster,SpatExtent-method}
\alias{relate,SpatRaster,SpatRaster-method}

\alias{is.related,SpatVector,SpatVector-method}
\alias{is.related,SpatVector,SpatExtent-method}
\alias{is.related,SpatExtent,SpatVector-method}
\alias{is.related,SpatVector,SpatRaster-method}
\alias{is.related,SpatRaster,SpatVector-method}
\alias{is.related,SpatRaster,SpatRaster-method}
\alias{is.related,SpatExtent,SpatRaster-method}
\alias{is.related,SpatRaster,SpatExtent-method}


\title{ 
Spatial relationships between geometries
}

\description{
\code{relate} returns a logical matrix indicating the presence or absence of a specific spatial relationships between the geometries in \code{x} and \code{y}. 

\code{is.related} returns a logical vector indicating the presence or absence of a specific spatial relationships between \code{x} and any of the geometries in \code{y}
}

\usage{
\S4method{relate}{SpatVector,SpatVector}(x, y, relation, pairs=FALSE, na.rm=TRUE)

\S4method{relate}{SpatVector,missing}(x, y, relation, pairs=FALSE, na.rm=TRUE)

\S4method{is.related}{SpatVector,SpatVector}(x, y, relation)
}

\arguments{
  \item{x}{SpatVector or SpatExtent}
  \item{y}{missing or as for \code{x}}
  \item{relation}{character. One of "intersects", "touches", "crosses", "overlaps", "within", "contains", "covers", "coveredby", "disjoint". Or a "DE-9IM" string such as "FF*FF****". See \href{https://en.wikipedia.org/wiki/DE-9IM}{wikipedia} or \href{https://docs.geotools.org/stable/userguide/library/jts/dim9.html}{geotools doc}}
  \item{pairs}{logical. If \code{TRUE} a two-column matrix is returned with the indices of the cases where the requested relation is \code{TRUE}. This is especially helpful when dealing with many geometries as the returned value is generally much smaller} 
  \item{na.rm}{logical. If \code{TRUE} and \code{pairs=TRUE}, geometries in \code{x} for which there is no related geometry in \code{y} are omitted}
} 

\value{
matrix (relate) or vector (is.related)
}

\seealso{
\code{\link{adjacent}}, \code{\link{nearby}}, \code{\link{intersect}}, \code{\link{crop}}
}

\examples{

# polygons
p1 <- vect("POLYGON ((0 0, 8 0, 8 9, 0 9, 0 0))")
p2 <- vect("POLYGON ((5 6, 15 6, 15 15, 5 15, 5 6))")
p3 <- vect("POLYGON ((8 2, 9 2, 9 3, 8 3, 8 2))")
p4 <- vect("POLYGON ((2 6, 3 6, 3 8, 2 8, 2 6))")
p5 <- vect("POLYGON ((2 12, 3 12, 3 13, 2 13, 2 12))")
p6 <- vect("POLYGON ((10 4, 12 4, 12 7, 11 7, 11 6, 10 6, 10 4))")

p <- rbind(p1, p2, p3, p4, p5, p6)
plot(p, col=rainbow(6, alpha=.5))
lines(p, lwd=2)
text(p)

## relate SpatVectors
relate(p1, p2, "intersects")
relate(p1, p3, "touches")
relate(p1, p5, "disjoint")
relate(rbind(p1, p2), p4, "disjoint")

## relate geometries within SpatVectors
# which are completely separated?
relate(p, relation="disjoint")

# which touch (not overlap or within)?
relate(p, relation="touches")
# which overlap (not merely touch, and not within)?
relate(p, relation="overlaps")
# which are within (not merely overlap)?
relate(p, relation="within")

# do they touch or overlap or are within?
relate(p, relation="intersects")

all(relate(p, relation="intersects") == 
  (relate(p, relation="overlaps") | 
   relate(p, relation="touches") |
   relate(p, relation="within")))

#for polygons, "coveredby" is "within"
relate(p, relation="coveredby")

# polygons, lines, and points 

pp <- rbind(p1, p2)
L1 <- vect("LINESTRING(1 11, 4 6, 10 6)")
L2 <- vect("LINESTRING(8 14, 12 10)")
L3 <- vect("LINESTRING(1 8, 12 14)")
lns <- rbind(L1, L2, L3)
pts <- vect(cbind(c(7,10,10), c(3,5,6)))

plot(pp, col=rainbow(2, alpha=.5))
text(pp, paste0("POL", 1:2), halo=TRUE)
lines(pp, lwd=2)
lines(lns, col=rainbow(3), lwd=4)
text(lns, paste0("L", 1:3), halo=TRUE)
points(pts, cex=1.5)
text(pts, paste0("PT", 1:3), halo=TRUE, pos=4)


relate(lns, relation="crosses")
relate(lns, pp, relation="crosses") 
relate(lns, pp, relation="touches") 
relate(lns, pp, relation="intersects") 

relate(lns, pp, relation="within")
# polygons can contain lines or points, not the other way around
relate(lns, pp, relation="contains") 
relate(pp, lns, relation="contains") 
# points and lines can be covered by polygons
relate(lns, pp, relation="coveredby")

relate(pts, pp, "within")
relate(pts, pp, "touches")
relate(pts, lns, "touches")

}


\keyword{methods}
\keyword{spatial}


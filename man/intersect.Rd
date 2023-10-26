\name{intersect}

\docType{methods}

\alias{intersect}

\alias{intersect,SpatVector,SpatVector-method}
\alias{intersect,SpatVector,SpatExtent-method}
\alias{intersect,SpatExtent,SpatVector-method}
\alias{intersect,SpatExtent,SpatExtent-method}
\alias{intersect,SpatRaster,SpatRaster-method}
\alias{intersect,SpatRaster,SpatExtent-method}
\alias{intersect,SpatExtent,SpatRaster-method}

\title{ 
Intersection
}

\description{
You can intersect SpatVectors with each other or with a SpatExtent. Intersecting points with points uses the extent of \code{y} to get the intersection. Intersecting of points and lines is not supported because of numerical inaccuracies with that. You can use \code{\link{buffer}}, to create polygons from lines and use these with intersect. 

You can also intersect two SpatExtents.

When intersecting two SpatRasters these need to be aligned (have the same origin and spatial resolution). The values of the returned SpatRaster are \code{TRUE} where both input rasters have values, \code{FALSE} where one has values, and \code{NA} in all other cells. 

When intersecting a SpatExtent and a SpatRaster, the SpatExtent is first aligned to the raster cell boundaries.

See \code{\link{crop}} for the intersection of a SpatRaster with a SpatExtent (or the extent of a SpatRaster or SpatVector) if you want a SpatRaster (not a SpatExtent) as output.

See \code{\link{is.related}(x, y, "intersects")} to find out which geometries of a SpatVector intersect. You can spatially subset a SpatVector with another one with \code{x\link{[}y]}.
}

\usage{
\S4method{intersect}{SpatVector,SpatVector}(x, y)

\S4method{intersect}{SpatVector,SpatExtent}(x, y)

\S4method{intersect}{SpatExtent,SpatVector}(x, y)

\S4method{intersect}{SpatExtent,SpatExtent}(x, y)

\S4method{intersect}{SpatRaster,SpatRaster}(x, y)

\S4method{intersect}{SpatRaster,SpatExtent}(x, y)

\S4method{intersect}{SpatExtent,SpatRaster}(x, y)
}

\arguments{
  \item{x}{SpatVector, SpatExtent, or SpatRaster}
  \item{y}{SpatVector, SpatExtent, or SpatRaster}
} 

\value{
Same as \code{x}
}

\seealso{
\code{\link{union}}, \code{\link{crop}}, \code{\link{relate}}, \code{\link{[}}
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

r <- s <- rast(ncol=5, nrow=5, xmin=1, xmax=5, ymin=1, ymax=5)
r[5:20] <- 5:20
s[11:20] <- 11:20
rs <- intersect(r, s)

u <- shift(r, .8)
us <- intersect(u, s)
}


\keyword{methods}
\keyword{spatial}


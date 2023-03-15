\name{sel}

\docType{methods}

\alias{sel}
\alias{sel,SpatRaster-method}
\alias{sel,SpatVector-method}

\title{
Spatial selection
}

\description{
Geometrically subset SpatRaster or SpatVector (to be done) by drawing on a plot (map).
}


\usage{
\S4method{sel}{SpatRaster}(x, ...)

\S4method{sel}{SpatVector}(x, use="rec", draw=TRUE, col="cyan", show=TRUE, ...)
}

\arguments{
\item{x}{SpatRaster or SpatVector}
\item{use}{character indicating what to draw. One of "rec" (rectangle) or "pol" (polygon)}
\item{draw}{logical. If \code{TRUE} the selected geometries are drawn on the map}
\item{col}{color to be used for drawing if \code{draw=TRUE}}
\item{show}{logical. If \code{TRUE} the area used to select geometries is drawn on the map}
\item{...}{additional graphics arguments for drawing the selected geometries}
}

\seealso{
\code{\link{crop}} and \code{\link{intersect}} to make an intersection and \code{\link{click}} and \code{\link{text}} to see cell values or geometry attributes.

Use \code{\link{draw}} to draw a SpatExtent of SpatVector that you want to keep.
}


\value{
SpatRaster or SpatVector
}


\examples{
\dontrun{
# select a subset of a SpatRaster
r <- rast(nrows=10, ncols=10)
values(r) <- 1:ncell(r)
plot(r)
s <- sel(r) # now click on the map twice

# plot the selection on a new canvas:
x11()
plot(s)

# vector
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
plot(v)
x <- sel(v) # now click on the map twice
x
}
}

\keyword{spatial}


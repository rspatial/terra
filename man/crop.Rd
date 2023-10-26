\name{crop}

\docType{methods}

\alias{crop}
\alias{crop,SpatRaster-method}
\alias{crop,SpatRasterDataset-method}
\alias{crop,SpatRasterCollection-method}
\alias{crop,SpatVector-method}
\alias{crop,SpatGraticule-method}


\title{Cut out a geographic subset}

\description{
Cut out a part of a SpatRaster or SpatVector. 

You can crop a SpatRaster with a SpatExtent, or with another object from which an extent can be obtained. Note that the SpatRaster returned may not have the exactly the same extent as the SpatExtent supplied because you can only select entire cells (rows and columns), and you cannot add new areas. See methods like \code{\link{resample}} and \code{\link{disagg}} to force SpatRasters to align and \code{\link{extend}} to add rows and/or columns.

You can only crop rectangular areas of a SpatRaster, but see argument \code{mask=TRUE} for setting cell values within SpatRaster to \code{NA}; or use the \code{\link{mask}} method after crop for additional masking options.

You can crop a SpatVector with another SpatVector. If these are not polygons, the minimum convex hull is used. Unlike with \code{\link{intersect}} the geometries and attributes of \code{y} are not transferred to the output. You can also crop a SpatVector with a rectangle (SpatRaster, SpatExtent).
}

\usage{
\S4method{crop}{SpatRaster}(x, y, snap="near", mask=FALSE, touches=TRUE, extend=FALSE, filename="", ...)

\S4method{crop}{SpatRasterDataset}(x, y, snap="near", extend=FALSE)

\S4method{crop}{SpatRasterCollection}(x, y, snap="near", extend=FALSE)

\S4method{crop}{SpatVector}(x, y, ext=FALSE)

\S4method{crop}{SpatGraticule}(x, y)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{y}{SpatRaster, SpatVector, SpatExtent, or any other object that has a SpatExtent (\code{\link{ext}} returns a \code{SpatExtent})}
  \item{snap}{character. One of "near", "in", or "out". Used to align \code{y} to the geometry of \code{x}} 
  \item{mask}{logical. Should \code{y} be used to mask? Only used if \code{y} is a SpatVector, SpatRaster or sf}  
  \item{touches}{logical. If \code{TRUE} and \code{mask=TRUE}, all cells touched by lines or polygons will be masked, not just those on the line render path, or whose center point is within the polygon}
  \item{extend}{logical. Should rows and/or columns be added if \code{y} is beyond the extent of \code{x}? Also see \code{\link{extend}}}  
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
  \item{ext}{logical. Use the extent of \code{y} instead of \code{y}. This also changes the behavior when \code{y} is an extent in two ways: (1) points that are on the extent boundary are removed and (2) lon/lat extents that go beyond -180 or 180 degrees longitude are wrapped around the earth to include areas at the other end of the dateline}
}


\value{
SpatRaster
}

\seealso{
\code{\link{intersect}}, \code{\link{extend}}

See \code{window} for a virtual and sometimes more efficient way to crop a dataset.
}

\examples{
r <- rast(xmin=0, xmax=10, ymin=0, ymax=10, nrows=25, ncols=25)
values(r) <- 1:ncell(r)
e <- ext(-5, 5, -5, 5)
rc <- crop(r, e)

# crop and mask
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
cm <- crop(r, v[9:12,], mask=TRUE)
plot(cm)
lines(v)

# crop vector
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
e <- ext(6.15, 6.3, 49.7, 49.8)
x <- crop(v, e)
plot(x, "NAME_1")
}

\keyword{spatial}


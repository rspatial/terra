\name{extractAlong}

\alias{extractAlong} 

\title{extract values along lines}

\description{ 
Extract raster values along a line. That is, the returned values are ordered along the line. That is not the case with \code{\link{extract}}
}

\usage{
extractAlong(x, y, ID=TRUE, cells=FALSE, xy=FALSE, online=FALSE, bilinear=TRUE)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatVector with lines geometry}
  \item{ID}{logical. Should an ID column be added? If so, the first column returned has the IDs (record numbers) of input SpatVector \code{y}}
  \item{cells}{logical. If \code{TRUE} the cell numbers are also returned}
  \item{xy}{logical. If \code{TRUE} the coordinates of the cells traversed by \code{y} are also returned. See \code{\link{xyFromCell}}}
  \item{online}{logical. If \code{TRUE} the returned coordinates are snapped to \code{y}} 
  \item{bilinear}{logical. If \code{TRUE} the returned raster values computed with bilinear interpolation from the nearest four cells. Only relevant if \code{online=TRUE}} 
}

\value{
data.frame
}

\seealso{
\code{\link{extract}}
}


\examples{
r <- rast(ncols=36, nrows=18, vals=1:(18*36))
cds1 <- rbind(c(-50,0), c(0,60), c(40,5), c(15,-45), c(-10,-25))
cds2 <- rbind(c(80,20), c(140,60), c(160,0), c(140,-55))
lines <- vect(list(cds1, cds2), "lines")

extractAlong(r, lines)
}

\keyword{methods}
\keyword{spatial}

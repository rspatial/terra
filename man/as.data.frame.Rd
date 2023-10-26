\name{as.data.frame}

\alias{as.data.frame}
\alias{as.data.frame,SpatRaster-method}
\alias{as.data.frame,SpatVector-method}
  

\title{SpatRaster or SpatVector to data.frame}

\description{
Coerce a SpatRaster or SpatVector to a data.frame
}

\usage{
\S4method{as.data.frame}{SpatVector}(x, row.names=NULL, optional=FALSE, geom=NULL, ...)

\S4method{as.data.frame}{SpatRaster}(x, row.names=NULL, optional=FALSE, xy=FALSE, cells=FALSE, na.rm=NA, ...)
}

\arguments{
 \item{x}{SpatRaster or SpatVector}
 \item{geom}{character or NULL. If not NULL, either "WKT" or "HEX", to get the geometry included in Well-Known-Text or hexadecimal notation. If \code{x} has point geometry, it can also be "XY" to add the coordinates of each point}
 \item{xy}{logical. If \code{TRUE}, the coordinates of each raster cell are included}
 \item{cells}{logical. If \code{TRUE}, the cell numbers of each raster cell are included}
 \item{na.rm}{logical. If \code{TRUE}, cells that have a NA value in at least one layer are removed. If the argument is set to \code{NA} only cells that have NA values in all layers are removed}
 \item{...}{Additional arguments passed to the \code{\link{data.frame}}}
 \item{row.names}{This argument is ignored}
 \item{optional}{This argument is ignored}
}
  
\seealso{\code{\link{as.list}, \link{as.matrix}}. See \code{\link{geom}} to only extract the geometry of a SpatVector}  
  
\value{
data.frame
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
as.data.frame(v)
}

\keyword{spatial}
\keyword{methods}


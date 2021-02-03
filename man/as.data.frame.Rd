\name{as.data.frame}

\alias{as.data.frame}
\alias{as.list}
\alias{as.data.frame,SpatRaster-method}
\alias{as.data.frame,SpatVector-method}
\alias{as.list,SpatVector-method}
  

\title{SpatRaster or SpatVector to data.frame}

\description{
Coerce a SpatRaster or SpatVector into a data.frame or a SpatVector into list
}

\usage{
\S4method{as.data.frame}{SpatVector}(x, geom=FALSE)

\S4method{as.data.frame}{SpatRaster}(x, xy=FALSE, cells=FALSE, na.rm=TRUE)

\S4method{as.list}{SpatVector}(x, geom=FALSE)
}

\arguments{
 \item{x}{SpatVector}
 \item{geom}{logical. If \code{TRUE} the WKT geometry is included}
 \item{xy}{logical. If \code{TRUE}, the coordinates of each raster cell are included}
 \item{cells}{logical. If \code{TRUE}, the cell numbers of each raster cell are included}
 \item{na.rm}{logical. If \code{TRUE}, cells that have a NA value in at least one layer are removed}
}
  
\seealso{see \code{\link{coerce}} for \code{as.data.frame} with a SpatRaster; and \code{\link{geom}} to only extract the geometry of a SpatVector}  
  
\value{
data.frame
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
as.data.frame(v)
as.list(v)
}

\keyword{spatial}
\keyword{methods}


\name{as.data.frame}

\alias{as.data.frame}
\alias{as.list}
\alias{as.data.frame,SpatRaster-method}
\alias{as.data.frame,SpatVector-method}
\alias{as.list,SpatRaster-method}
\alias{as.list,SpatVector-method}
  

\title{SpatRaster or SpatVector to data.frame}

\description{
Coerce a SpatRaster or SpatVector to a data.frame or coerce a SpatVector to a list
}

\usage{
\S4method{as.data.frame}{SpatVector}(x, row.names=NULL, optional=FALSE, geom=NULL, ...)

\S4method{as.data.frame}{SpatRaster}(x, row.names=NULL, optional=FALSE, xy=FALSE, cells=FALSE, na.rm=TRUE, ...)

\S4method{as.list}{SpatVector}(x, geom=NULL)

\S4method{as.list}{SpatRaster}(x)
}

\arguments{
 \item{x}{SpatRaster or SpatVector}
 \item{geom}{character or NULL. If not NULL, either "WKT" or "HEX", to get the geometry included in Well-Known-Text or hexadecimal notation}
 \item{xy}{logical. If \code{TRUE}, the coordinates of each raster cell are included}
 \item{cells}{logical. If \code{TRUE}, the cell numbers of each raster cell are included}
 \item{na.rm}{logical. If \code{TRUE}, cells that have a NA value in at least one layer are removed}
 \item{...}{Additional arguments passed to the \code{\link{data.frame}}}
 \item{row.names}{This argument is ignored}
 \item{optional}{This argument is ignored}
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


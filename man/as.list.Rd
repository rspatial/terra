\name{as.list}

\alias{as.list}
\alias{as.list,SpatRaster-method}
\alias{as.list,SpatRasterDataset-method}
\alias{as.list,SpatRasterCollection-method}
\alias{as.list,SpatVector-method}
\alias{as.list,SpatVectorCollection-method}
  

\title{Coerce a Spat* object to a list}

\description{
Coerce a SpatRaster, SpatRasterCollection, SpatRasterDataset, SpatVector or SpatVectorCollection to a list. With a SpatRaster, each layer becomes a list element. With a SpatRasterCollection or SpatRasterDataset, each SpatRaster becomes a list element. With a SpatVector, each variable (attribute) becomes a list element. With a SpatVectorCollection, each SpatVector becomes a list element.
}

\usage{
\S4method{as.list}{SpatRaster}(x, ...)

\S4method{as.list}{SpatRasterCollection}(x, ...)

\S4method{as.list}{SpatVector}(x, geom=NULL, ...)

\S4method{as.list}{SpatVectorCollection}(x, ...)
}

\arguments{
 \item{x}{SpatRaster, SpatRasterDataset, SpatRasterCollection, or SpatVector}
 \item{geom}{character or NULL. If not NULL, either "WKT" or "HEX", to get the geometry included in Well-Known-Text or hexadecimal notation. If \code{x} has point geometry, it can also bey "XY" to add the coordinates of each point}
 \item{...}{Additional arguments. These are ignored}
}
  
\seealso{see \code{\link{coerce}} for \code{as.data.frame} with a SpatRaster; and \code{\link{geom}} to only extract the geometry of a SpatVector}  
  
\value{
list
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
as.list(v)


s <- rast(system.file("ex/logo.tif", package="terra")) + 1  
as.list(s)
}

\keyword{spatial}
\keyword{methods}


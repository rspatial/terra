\name{wrap}

\alias{wrap}
\alias{unwrap}

\alias{wrap,SpatVector-method}
\alias{wrap,SpatRaster-method}
\alias{wrap,SpatRasterDataset-method}
\alias{wrap,SpatRasterCollection-method}
\alias{unwrap,ANY-method}
\alias{unwrap,PackedSpatRaster-method}
\alias{unwrap,PackedSpatRasterDC-method}
\alias{unwrap,PackedSpatVector-method}


\title{wrap and unwrap}

\description{
Use \code{wrap} to pack a SpatVector or SpatRaster* to create a Packed* object. Packed objects can be passed over a connection that serializes (e.g. to nodes on a computer cluster). At the receiving end they need to be unpacked with \code{unwrap}.
}

\usage{
\S4method{wrap}{SpatRaster}(x, proxy=FALSE, path=NULL, overwrite=FALSE)

\S4method{wrap}{SpatRasterDataset}(x, proxy=FALSE)

\S4method{wrap}{SpatRasterCollection}(x, proxy=FALSE)

\S4method{wrap}{SpatVector}(x)

\S4method{unwrap}{ANY}(x)
}

\arguments{
  \item{x}{SpatVector, SpatRaster, SpatRasterDataset or SpatRasterCollection}
  \item{proxy}{logical. If \code{FALSE} raster cell values are forced to memory if possible. If \code{TRUE}, a reference to source filenames is stored for data sources that are not in memory}
  \item{path}{character. If not \code{NULL}, the path where raster files are copied to. This should normally not be used}
  \item{overwrite}{Should existing files be overwritten when \code{path} is not \code{NULL}? If this value is not TRUE or FALSE, only files that do not exist are copied}
}

\value{
\code{wrap}: Packed* object

\code{unwrap}: SpatVector, SpatRaster, SpatRasterCollection, SpatRasterDataset
}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
p <- wrap(v)
p
vv <- vect(p)
vv
}


\keyword{ spatial }
\keyword{ methods }


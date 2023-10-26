\name{serialize}

\alias{serialize}
\alias{saveRDS}
\alias{unserialize}
\alias{readRDS}

\alias{serialize,SpatVector-method}
\alias{serialize,SpatRaster-method}
\alias{serialize,SpatRasterDataset-method}
\alias{serialize,SpatRasterCollection-method}
\alias{unserialize,ANY-method}

\alias{saveRDS,SpatVector-method}
\alias{saveRDS,SpatRaster-method}
\alias{saveRDS,SpatRasterCollection-method}
\alias{saveRDS,SpatRasterDataset-method}
\alias{readRDS,character-method}


\title{saveRDS and serialize for SpatVector and SpatRaster*}

\description{
serialize and saveRDS for SpatVector, SpatRaster, SpatRasterDataset and SpatRasterCollection. Note that these objects will first be "packed" with \code{\link{wrap}}, and after unserialize/readRDS they need to be unpacked with \code{rast} or \code{vect}.

Extensive use of these functions is not recommended. Especially for SpatRaster it is generally much more efficient to use \code{\link{writeRaster}} and write, e.g., a GTiff file. 
}

\usage{
\S4method{saveRDS}{SpatRaster}(object, file="", ascii = FALSE, version = NULL, compress=TRUE, refhook = NULL)

\S4method{saveRDS}{SpatRasterDataset}(object, file="", ascii = FALSE, version = NULL, compress=TRUE, refhook = NULL)

\S4method{saveRDS}{SpatRasterCollection}(object, file="", ascii = FALSE, version = NULL, compress=TRUE, refhook = NULL)

\S4method{saveRDS}{SpatVector}(object, file="", ascii = FALSE, version = NULL, compress=TRUE, refhook = NULL)

\S4method{serialize}{SpatRaster}(object, connection, ascii = FALSE, xdr = TRUE, version = NULL, refhook = NULL)

\S4method{serialize}{SpatVector}(object, connection, ascii = FALSE, xdr = TRUE, version = NULL, refhook = NULL)
}

\arguments{
  \item{object}{SpatVector, SpatRaster, SpatRasterDataset or SpatRasterCollection}
  \item{file}{file name to save object to}
  \item{connection}{see \code{\link{serialize}}}
  \item{ascii}{see \code{\link{serialize}} or \code{\link{saveRDS}}}
  \item{version}{see \code{\link{serialize}} or \code{\link{saveRDS}}}
  \item{compress}{see \code{\link{serialize}} or \code{\link{saveRDS}}}
  \item{refhook}{see \code{\link{serialize}} or \code{\link{saveRDS}}}
  \item{xdr}{see \code{\link{serialize}} or \code{\link{saveRDS}}}
}

\value{
Packed* object 
}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
p <- serialize(v, NULL)
head(p)
x <- unserialize(p)
x
}


\keyword{ spatial }
\keyword{ methods }


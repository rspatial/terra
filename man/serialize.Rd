\name{serialize}

\alias{serialize}
\alias{saveRDS}
\alias{unserialize}
\alias{readRDS}

\alias{serialize,SpatVector-method}
\alias{serialize,SpatRaster-method}
\alias{unserialize,ANY-method}

\alias{saveRDS,SpatVector-method}
\alias{saveRDS,SpatRaster-method}
\alias{readRDS,character-method}


\title{serialize and saveRDS for SpatRaster and SpatVector}

\description{
serialize and saveRDS for SpatRaster and SpatVector. Note that these objects will first be "packed" with \code{\link{wrap}}, and after unserialize/readRDS they need to be unpacked with \code{rast} or \code{vect}.

Use of these functions is not recommended. Especially for SpatRaster it is generally much more efficient to use \code{\link{writeRaster}} and write, e.g., a GTiff file. 

SpatRaster objects must have all values in memory (that is, the cell values are not in files) to be serialized. These functions use \code{\link{set.values}} to load values into memory if needed and if deemed possible given the amount of RAM available. 
}

\usage{
\S4method{saveRDS}{SpatRaster}(object, file="", ascii = FALSE, version = NULL, compress=TRUE, refhook = NULL)

\S4method{saveRDS}{SpatVector}(object, file="", ascii = FALSE, version = NULL, compress=TRUE, refhook = NULL)

\S4method{serialize}{SpatRaster}(object, connection, ascii = FALSE, xdr = TRUE, version = NULL, refhook = NULL)

\S4method{serialize}{SpatVector}(object, connection, ascii = FALSE, xdr = TRUE, version = NULL, refhook = NULL)
}

\arguments{
  \item{object}{SpatVector or SpatRaster}
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


\name{serialize}

\alias{serialize}
\alias{saveRDS}

\alias{serialize,SpatVector-method}
\alias{serialize,SpatRaster-method}
\alias{saveRDS,SpatVector-method}
\alias{saveRDS,SpatRaster-method}


\title{serialize and saveRDS for SpatRaster and SpatVector}

\description{
serialize or saveRDS for SpatRaster and SpatVector. Note that these objects will first be "packed" with \code{\link{wrap}}, and after unserialize/readRDS they need to be unpacked with \code{rast} or \code{vect}.

Only SpatRaster objects that have all values in memory can be serialized. That is, SpatRaster that have (some of) their values in a file cannot be serialized. 
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
vect(x)
}


\keyword{ spatial }
\keyword{ methods }


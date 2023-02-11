\name{wrap}

\alias{wrap}
\alias{unwrap}

\alias{wrap,SpatVector-method}
\alias{wrap,SpatRaster-method}
\alias{unwrap,ANY-method}
\alias{unwrap,PackedSpatRaster-method}
\alias{unwrap,PackedSpatVector-method}


\title{wrap and unwrap}

\description{
Use \code{wrap} to pack a SpatRaster or SpatVector object to create a Packed* object. Packed objects can be passed over a connection that serializes (e.g. to nodes on a computer cluster). At the receiving end they need to be unpacked with \code{unwrap}.
}

\usage{
\S4method{wrap}{SpatRaster}(x, proxy=FALSE)

\S4method{wrap}{SpatVector}(x)

\S4method{unwrap}{ANY}(x)
}

\arguments{
  \item{x}{SpatVector or SpatRaster}
  \item{proxy}{logical. If \code{FALSE} raster cell values are forced to memory if possible. If \code{TRUE}, a reference to source filenames is stored for data sources that are not in memory}
}

\value{
\code{wrap}: Packed* object

\code{unwrap}: SpatRaster or SpatVector
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


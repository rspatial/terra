\name{wrap}

\alias{wrap}

\alias{wrap,SpatVector-method}
\alias{wrap,SpatRaster-method}


\title{wrap (pack) a SpatRaster or SpatVector object}

\description{
Wrap a SpatRaster or SpatVector object to create a Packed* object. Packed objects can be passed over a connection that serializes (e.g. to nodes on a computer cluster). At the receiving end they need to be unpacked with either \code{rast} or \code{vect}.
}

\usage{
\S4method{wrap}{SpatRaster}(x, proxy=FALSE)

\S4method{wrap}{SpatVector}(x)
}

\arguments{
  \item{x}{SpatVector or SpatRaster}
  \item{proxy}{logical. If \code{FALSE} raster cell values are forced to memory if possible. If \code{TRUE}, a reference to source filenames is stored for data sources that are not in memory}
}

\value{
Packed* object 
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


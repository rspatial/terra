\name{wrap}

\alias{wrap}

\alias{wrap,SpatVector-method}
\alias{wrap,SpatRaster-method}
\alias{wrap,Spatial-method}


\title{wrap (pack) a SpatRaster or SpatVector object}

\description{
Wrap a SpatRaster or SpatVector object to create a Packed* object. Packed objects can be saved as an R object to disk (.rds or .RData), or passed over a connection that serializes (e.g. to nodes on a computer cluster); but with large datasets passing a filename could be more sensible in that context. 
}

\usage{
\S4method{wrap}{SpatRaster}(x)

\S4method{wrap}{SpatVector}(x)
}

\arguments{
  \item{x}{SpatVector or SpatRaster}
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


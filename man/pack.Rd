\name{wrap}

\alias{wrap,SpatVector-method}
\alias{wrap,SpatRaster-method}
\alias{wrap,Spatial-method}
\alias{pack,missing-method}

\alias{pack}

\title{pack a SpatRaster or SpatVector object}

\description{
Wrap a SpatRaster or SpatVector object to create a Packed* object that it can be saved as an R object to disk, or passed over a connection that serializes (e.g. using a computer cluster). 
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


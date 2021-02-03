\name{pack}

\alias{pack,SpatVector-method}
\alias{pack,SpatRaster-method}
\alias{pack,Spatial-method}

\alias{pack}

\title{pack a SpatRaster or SpatVector object}

\description{
Pack a SpatRaster or SpatVector object so that it can be saved as an R object to disk, or passed over a connection that serializes (e.g. using a computer cluster). 
}

\usage{
\S4method{pack}{SpatRaster}(x)

\S4method{pack}{SpatVector}(x)
}

\arguments{
  \item{x}{SpatVector or SpatRaster}
  \item{...}{additional arguments. None implemented}
}

\value{
Packed* object 
}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
p <- pack(v)
p
vv <- vect(p)
vv
}


\keyword{ spatial }
\keyword{ methods }


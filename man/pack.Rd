\name{pack}

\alias{pack,SpatVector-method}
\alias{pack,SpatRaster-method}
\alias{pack,Spatial-method}

\alias{pack}

\title{pack a Spat* object}

\description{
Pack a Spat* object so that it can be saved as an R object to disk, or passed over a connection that serializes (e.g. using a computer cluster). 
}

\usage{
\S4method{pack}{SpatRaster}(x, ...)

\S4method{pack}{SpatVector}(x, ...)
}

\arguments{
  \item{x}{SpatVector or SpatRaster}
  \item{...}{additional arguments. None implemented}
}

\value{
Packed* object 
}


\keyword{ spatial }
\keyword{ methods }


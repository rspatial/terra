\name{as.data.frame}

\alias{as.data.frame}
\alias{as.data.frame,SpatLayer-method}
  

\title{Get the attributes of a  SpatLayer}

\description{
Get a data.frame with the attribute values of a SpatLayer
}

\usage{
\S4method{as.data.frame}{SpatLayer}(x, row.names=NULL, optional=FALSE, ...)
}

\arguments{
 \item{x}{SpatLayer}
 \item{row.names}{ignored} 
 \item{optional}{ignored}
 \item{...}{Additional arguments (none)} 
}
  
\value{
data.frame
}

\examples{
f <- system.file("exdata/lux.shp", package="terra")
v <- vect(f)
as.data.frame(v)
}

\keyword{spatial}
\keyword{methods}

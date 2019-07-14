\name{as.data.frame}

\alias{as.data.frame}
\alias{as.data.frame,SpatVector-method}
  

\title{Get the attributes of a  SpatVector}

\description{
Get a data.frame with the attribute values of a SpatVector
}

\usage{
\S4method{as.data.frame}{SpatVector}(x, row.names=NULL, optional=FALSE, ...)
}

\arguments{
 \item{x}{SpatVector}
 \item{row.names}{ignored} 
 \item{optional}{ignored}
 \item{...}{additional arguments (none)} 
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

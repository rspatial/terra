\name{add}

\docType{methods}

\alias{add<-}
\alias{add<-,SpatRaster,SpatRaster-method}

\title{Add (in place) a SpatRaster to another SpatRaster object}

\description{
Add (in place) a SpatRaster to another SpatRaster object. Comparable with \code{\link{c}}, but withouth copying the object. 
}

\usage{
\S4method{add}{SpatRaster,SpatRaster}(x)<-value
}

\arguments{
  \item{x}{SpatRaster}
  \item{value}{SpatRaster}
}

\seealso{\code{\link{c}}}

\value{
SpatRaster
}

\examples{
r <- rast(nrow=5, ncol=9, vals=1:45)
x <- c(r, r*2)
add(x) <- r*3
x
}

\keyword{spatial}

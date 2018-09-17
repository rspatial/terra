\name{values}

\docType{methods}


\alias{values}
\alias{values,SpatRaster-method}
\alias{values<-}
\alias{values<-,SpatRaster,numeric-method}
  
\title{Get cell values}

\description{
values returns all cell values of a SpatRaster object.
}

\usage{
\S4method{values}{SpatRaster}(x,...)
\S4method{values}{SpatRaster,numeric}(x,...)<-value
}

\arguments{
 \item{x}{SpatRaster object}
 \item{value}{numeric}
 \item{...}{additional arguments. none implemented}
}

\value{
matrix
}


\examples{
r <- rast(system.file("external/test.grd", package="terra"))
r
v <- values(r)
values(r) <- v * 10
r
}

\keyword{spatial}
\keyword{methods}

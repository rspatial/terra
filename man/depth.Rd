\name{depth}

\alias{depth}
\alias{depth<-}
\alias{depth,SpatRaster-method}
\alias{depth<-,SpatRaster-method}
\alias{depthName}
\alias{depthName<-}
\alias{depthName,SpatRaster-method}
\alias{depthName<-,SpatRaster-method}

  
\title{depth of SpatRaster layers}

\description{
Get or set the depth of the layers of a SpatRaster. Experimental. 
}

\usage{
\S4method{depth}{SpatRaster}(x)

\S4method{depth}{SpatRaster}(x)<-value

\S4method{depthName}{SpatRaster}(x)

\S4method{depthName}{SpatRaster}(x)<-value
}

\arguments{
  \item{x}{SpatRaster}
  \item{value}{numeric vector (\code{depth}), or character (\code{depthName})}
}

\value{
numeric  
}

\seealso{\code{\link{time}}}

\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   

depth(s) <- 1:3
depth(s)
}

\keyword{spatial}


\name{plot_extent}

\docType{methods}

\alias{plot,SpatExtent,missing-method}

\title{Plot a SpatExtent}

\description{
Plot a SpatExtent. Use \code{\link{lines}} to add a SpatExtent to an existing map.

See \code{\link{plot}} for plotting other object types.
}

\usage{
\S4method{plot}{SpatExtent,missing}(x, y, ...)
}

\arguments{
  \item{x}{SpatExtent}
  \item{y}{missing}
  \item{...}{additional graphical arguments for lines}  
}

\seealso{
\code{\link{plot}}
}	

\examples{ 
r <- rast()
plot(ext(r))
}


\keyword{methods}
\keyword{spatial}


\name{plot_extent}

\docType{methods}

\alias{plot,SpatExtent,missing-method}

\title{Plot a SpatExtent}

\description{
Plot a SpatExtent. Use \code{\link{lines}} to add a SpatExtent to an existing map.

There are other help files for plotting a \code{\link[terra:plot]{SpatRaster}}, \code{\link[terra:plot_vector]{SpatVector}} or \code{\link[terra:plot_graticule]{SpatGraticule}}
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
plot a \code{\link[terra:plot]{SpatRaster}}, \code{\link[terra:plot_vector]{SpatVector}}, or a \code{\link[terra:plot_graticule]{SpatGraticule}}.

\code{\link{points}}, \code{\link{lines}}, \code{\link{polys}}, 
scale bar: \code{\link{sbar}}, north arrow: \code{\link{north}}
}	

\examples{ 
r <- rast()
plot(ext(r))
}


\keyword{methods}
\keyword{spatial}


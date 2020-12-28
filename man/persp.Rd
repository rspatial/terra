\name{persp}

\docType{methods}

\alias{persp}
\alias{persp,SpatRaster-method}

\title{Perspective plot}

\description{
Perspective plot of a SpatRaster. This is an implementation of a generic function in the graphics package.
}


\usage{
\S4method{persp}{SpatRaster}(x, maxcells=100000, ...)
}

\arguments{
  \item{x}{SpatRaster. Only the first layer is used}
  \item{maxcells}{integer > 0. Maximum number of cells to use for the plot. If \code{maxpixels < ncell(x)}, \code{spatSample(method="regular")} is used before plotting}
  \item{...}{Any argument that can be passed to \code{\link[graphics]{persp}} (graphics package)}  
}


\seealso{ \code{\link[graphics]{persp}}, \code{contour}, \code{plot} }



\examples{ 
r <- rast(system.file("ex/elev.tif", package="terra"))
persp(r)
}

\keyword{methods}
\keyword{spatial}


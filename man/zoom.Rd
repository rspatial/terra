\name{zoom}

\docType{methods}

\alias{zoom}
\alias{zoom,SpatRaster-method}
\alias{zoom,SpatVector-method}


\title{Zoom in on a map}


\description{
Zoom in on a map (plot) by providing a new extent, by default this is done by clicking twice on the map.
}

\usage{
\S4method{zoom}{SpatRaster}(x, e=draw(), maxcell=100000, layer=1, new=FALSE, ...)

\S4method{zoom}{SpatVector}(x, e=draw(), new=FALSE, ...)
}

\arguments{
\item{x}{SpatRaster}
\item{e}{SpatExtent}
\item{maxcell}{positive integer. Maximum number of cells used for the map}
\item{layer}{positive integer to select the layer to be used}
\item{new}{logical. If \code{TRUE}, the zoomed in map will appear on a new device (window)}
\item{...}{additional arguments passed to \code{\link{plot}}}
}

\value{
SpatExtent (invisibly)
}


\seealso{ \code{\link{draw}}, \code{\link{plot}}}

\keyword{spatial}


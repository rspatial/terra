\name{zoom}

\docType{methods}

\alias{zoom}
\alias{zoom,SpatRaster-method}


\title{Zoom in on a map}


\description{
Zoom in on a map (plot) by providing a new extent, by default this is done by clicking twice on the map.
}

\usage{
\S4method{zoom}{SpatRaster}(x, ext=draw(), maxcell=10000, layer=1, new=TRUE, ...)
}

\arguments{
\item{x}{SpatRaster}
\item{ext}{SpatExtent}
\item{maxcell}{Maximum number of cells used for the map}
\item{layer}{Positive integer to select the layer to be used}
\item{new}{Logical. If \code{TRUE}, the zoomed in map will appear on a new device (window)}
\item{...}{additional paramters for \code{\link[graphics]{plot}}}
}

\value{
SpatExtent (invisibly)
}


\seealso{ \code{\link{draw}}, \code{\link{plot}}}

\keyword{spatial}


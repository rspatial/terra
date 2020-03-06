\name{select}

\docType{methods}

\alias{select}
\alias{select,SpatRaster-method}

\title{
Geometric subsetting
}

\description{
Geometrically subset SpatRaster or SpatVector (to be done) by drawing on a plot (map).
}


\usage{
\S4method{select}{SpatRaster}(x, ...)
}

\arguments{
\item{x}{SpatRaster or SpatVector}
\item{...}{additional arguments passed on to \code{\link{draw}}}
}

\seealso{
\code{\link{click}, \link{crop}}
}


\value{
SpatRaster or SpatVector
}


\examples{
\donttest{
# select a subset of a RasterLayer
r <- rast(nrow=10, ncol=10)
values(r) <- 1:ncell(r)
plot(r)
s <- select(r) # now click on the map twice

# plot the selection on a new canvas:
x11()
plot(s)
}
}

\keyword{spatial}


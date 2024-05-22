\name{add_box}

\alias{add_box}

\title{draw a box}

\description{ 
Similar to \code{\link[graphics]{box}} allowing adding a box around a map. This function will place the legend in the locations within the mapped area as delineated by the axes.
}

\usage{
add_box(...)
}

\arguments{
  \item{...}{arguments passed to \code{\link{lines}}}
}


\seealso{\code{\link{add_legend}}, \code{\link{add_grid}}}

\examples{
v <- vect(system.file("ex/lux.shp", package="terra"))
plot(v)
add_box(col="red", lwd=3, xpd=TRUE)
}

\keyword{methods}
\keyword{spatial}


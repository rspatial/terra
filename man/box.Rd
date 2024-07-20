\name{add_box}

\alias{add_box}

\title{draw a box}

\description{ 
Similar to \code{\link[graphics]{box}} allowing adding a box around a map. This function will place the box around the mapped area.
}

\usage{
add_box(...)
}

\arguments{
  \item{...}{arguments passed to \code{\link{lines}}}
}


\seealso{\code{\link{add_legend}}, \code{\link{add_grid}}, \code{\link{add_mtext}}}

\examples{
v <- vect(system.file("ex/lux.shp", package="terra"))
plot(v)
add_box(col="red", lwd=3, xpd=TRUE)
}

\keyword{methods}
\keyword{spatial}


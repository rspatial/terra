\name{add_legend}

\alias{add_legend}

\title{add a custom legend}

\description{ 
Wrapper around \code{\link{legend}} that allows adding a custom legend to a map using a keyword such as "topleft" or "bottomright". This function will place the legend in the locations within the mapped area as delineated by the axes.
}

\usage{
add_legend(x, y, ...)
}


\arguments{
  \item{x}{The keyword to be used to position the legend (or the x coordinate)}
  \item{y}{The y coordinate to be used to position the legend (is x is also a coordinate)}
  \item{...}{arguments passed to \code{\link{legend}}}
}


\seealso{\code{\link{add_box}}, \code{\link{add_grid}}}

\examples{
v <- vect(system.file("ex/lux.shp", package="terra"))
plot(v)
points(centroids(v), col="red")
legend("topleft", legend = "centroids", pch = 20, xpd=NA, bg="white", col="red")
add_legend("topright", legend = "centroids", pch = 20, col="red")
}

\keyword{methods}
\keyword{spatial}


\name{add_grid}

\alias{add_grid}

\title{add a grid to a map made with terra}

\description{ 
Adaptation of \code{\link[graphics]{grid}} that allows adding a grid to a map. This function will place the legend in the locations within the mapped area as delineated by the axes.

Also see \code{\link{graticule}}
}

\usage{
add_grid(nx=NULL, ny=nx, col="lightgray", lty="dotted", lwd=1)
}

\arguments{
  \item{nx, ny}{number of cells of the grid in x and y direction. When NULL, as per default, the grid aligns with the tick marks on the corresponding default axis (i.e., tickmarks as computed by axTicks). When NA, no grid lines are drawn in the corresponding direction}
  \item{col}{character or (integer) numeric; color of the grid lines}
  \item{lty}{character or (integer) numeric; line type of the grid lines}
  \item{lwd}{non-negative numeric giving line width of the grid lines}
}

\seealso{ \code{\link{graticule}}, \code{\link{add_box}}, \code{\link{add_grid}}}

\examples{
v <- vect(system.file("ex/lux.shp", package="terra"))
plot(v)
add_grid()
}

\keyword{methods}
\keyword{spatial}


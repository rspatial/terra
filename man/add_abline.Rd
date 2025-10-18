\name{add_abline}

\alias{add_abline}

\title{add vertical and/or horizontal lines to a map made with terra}

\description{ 
Adaptation of \code{\link[graphics]{abline}} that allows adding a horizonal or vertical lines to a map. This function will place the lines in the locations within the mapped area as delineated by the axes. It is meant to be used when you specifiy your own tick marks, such that \code{\link{add_grid}} does not work.

Also see \code{\link{graticule}}
}

\usage{
add_abline(h=NULL, v=NULL, ...)
}

\arguments{
  \item{h}{the y-value(s) for horizontal line(s)}
  \item{v}{the x-value(s) for vertical line(s)}
  \item{...}{additional graphical parameters for drawing lines}
}

\seealso{\code{\link{add_grid}}, \code{\link{graticule}}, \code{\link{add_legend}}, \code{\link{add_box}}, \code{\link{add_grid}}, \code{\link{add_mtext}}}

\examples{
v <- vect(system.file("ex/lux.shp", package="terra"))
atx <- seq(xmin(v), xmax(v), .1)
aty <- seq(ymin(v), ymax(v), .1)
plot(v, pax=list(xat=atx, yat=aty), ext=ext(v)+.2)
add_abline(h=aty, v=atx, lty=2, col="gray")
}

\keyword{methods}
\keyword{spatial}


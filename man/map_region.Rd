\name{map_region}

\alias{map_region}

\title{Get the coordinates of the extent of a map}

\description{ 
Helper function for creating custom map elements that are aligned with the axes of a map (base plot created with a SpatRaster and/or SpatVector). For example, you may need to know the coordinates for the upper-left corner of a map to add some information there. 

Unlike the standard base plot, terra keeps the axis aligned with the data. For that reason you cannot use \code{par()$usr} to get these coordinates.

The coordinates returned by this function are used in, for example, code{\link{add_legend}} such that a legend can be automatically placed in the a particular corner.

This function only returns meaningful results of the active plot (canvas) was create with a call to \code{plot} with a SpatRaster or SpatVector as first argument.
}

\usage{
map_region(...)
}

\seealso{\code{\link{add_legend}}, \code{\link{add_grid}}, \code{\link{add_box}}}

\examples{
r <- rast(xmin=0, xmax=10, ymin=0, ymax=10, res=1, vals=1:100)
plot(r)

map_region()
par()$usr
}

\keyword{methods}
\keyword{spatial}


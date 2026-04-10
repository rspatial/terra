\name{legend_cont}

\alias{legend_cont}

\title{Add a continuous legend}

\description{
Add a continuous (color gradient) legend to an existing plot. This can be used with any base plot, not just \code{\link{plot,SpatRaster,numeric-method}}. The legend is drawn as a color bar with tick marks and labels.
}

\usage{
legend_cont(x="right", y=NULL, legend, col,
    size=NULL, title=NULL, at=NULL, digits=NULL, ...)
}


\arguments{
  \item{x}{character or numeric. Position of the legend. Keywords: \code{"right"}, \code{"left"}, \code{"top"}, \code{"bottom"}, \code{"topright"}, \code{"bottomright"}. Or a numeric x-coordinate (see \code{y})}
  \item{y}{numeric. Optional y-coordinate(s) for the legend position, used when \code{x} is numeric}
  \item{legend}{SpatRaster or numeric vector. Used to determine the value range of the legend. If a SpatRaster, the range is taken from \code{\link{minmax}}. If numeric (e.g. \code{c(min, max)}), the range of the values is used}
  \item{col}{character. Vector of colors for the gradient. If missing, the default palette is used (\code{map.pal("viridis", 100)})}
  \item{size}{numeric. One or two values to control the size of the legend bar. The first value scales the length (0 to 1), the second scales the width}
  \item{title}{character. Title to display above or beside the legend bar}
  \item{at}{numeric. Specific values at which to place tick marks. If \code{NULL} (the default), ticks are placed automatically with \code{\link{pretty}}}
  \item{digits}{non-negative integer. Number of decimal places for tick labels. If \code{NULL} (the default) an appropriate number is computed from the range}
  \item{...}{Additional legend parameters: \code{bg} (background color behind the legend, e.g. \code{"white"}), \code{cex}, \code{horiz}, \code{reverse}, \code{labels}, \code{tic} (or \code{tick}), \code{tic.col}, \code{tic.lwd}, \code{tic.box.col}, and others (see the \code{plg} argument in \code{\link{plot,SpatRaster,numeric-method}} for details)}
}


\value{
Invisible list with the legend parameters (can be used for further customization).
}


\seealso{\code{\link{add_legend}}, \code{\link{add_box}}, \code{\link{plot,SpatRaster,numeric-method}}}


\examples{
# with a SpatRaster
r <- rast(ncols=40, nrows=40, xmin=0, xmax=1, ymin=0, ymax=1, vals=runif(1600))
plot(r, legend=FALSE, mar=c(3.1, 3.1, 2.1, 7.1))
legend_cont("right", legend=r)

# horizontal
plot(r, legend=FALSE, mar=c(6.1, 3.1, 2.1, 2.1))
legend_cont("bottom", legend=r, horiz=TRUE, title="values")

# on the map 
plot(r, legend=FALSE, alpha=0.1)
legend_cont(0.1, legend=r)
legend_cont(0.3, c(0.5, 0.9), legend=r, bg="white")
legend_cont(c(0.4, 0.8), 0.2, horiz=TRUE, legend=r, bg="white")


# with a numeric range and custom colors
cols <- heat.colors(100)
vals <- seq(5,95,10)
par(mar=c(3.1, 3.1, 2.1, 7.1))
plot(1:10, col=cols[vals], cex=2, pch=20)
legend_cont("right", legend=c(1, 100), col=cols)

}

\keyword{methods}
\keyword{spatial}


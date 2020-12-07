\name{plotVector}

\docType{methods}

\alias{plot,SpatVector,missing-method}
\alias{plot,SpatVector,character-method}
\alias{plot,SpatVector,numeric-method}
\alias{plot,SpatExtent,missing-method}


\title{Plot a SpatVector}

\description{
Plot SpatVector data to make a map
}

\usage{

\S4method{plot}{SpatVector,character}(x, y, col, type,
	mar=NULL, legend=TRUE, add=FALSE, axes=!add, main=y,
    plg=list(), pax=list(), nr, nc, ...)

\S4method{plot}{SpatVector,numeric}(x, y, ...)

\S4method{plot}{SpatVector,missing}(x, y, ...)

\S4method{plot}{SpatExtent,missing}(x, y, ...)
}

\arguments{
  \item{x}{SpatVector}
  \item{y}{missing or positive integer indicating the layer(s) to be plotted, or the name of the layer}
  \item{col}{character. Colors}
  \item{type}{character. Type of map/legend. One of "none", "continuous", "classes", "interval"} 
  \item{mar}{numeric vector of lenght 4 to set the margins of the plot (to make space for the legend). The default is c(3.1, 3.1, 2.1, 7.1) for a plot with a legend. Use \code{mar=NA} to not set the margins}
  \item{legend}{logical. Draw a legend?}
  \item{add}{logical. If \code{TRUE} add the object to the current plot}
  \item{axes}{logical. Draw axes?}
  \item{main}{character. Text for the main title of the plot}
  \item{plg}{list with parameters for drawing the legend. See the arguments for \code{\link{legend}}}
  \item{pax}{list with parameters for drawing axes. See the arguments for \code{\link{axis}}}  
  \item{nc}{positive integer. Optional. The number of columns to divide the plotting device in (when plotting multiple layers)}
  \item{nr}{positive integer. Optional. The number of rows to divide the plotting device in (when plotting multiple layers)}
  \item{...}{additional arguments passed to \code{plot("SpatVector", "character")} method, and graphical arguments such as \code{lwd}, \code{cex} and \code{pch}}  
}

\seealso{\code{\link{sbar}}}

\examples{ 
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)

plot(v)

plot(v, 2, pax=list(sides=1:2), plg=list(x=6.2, y=50.2, cex=1.2))

plot(v, 4, pax=list(sides=1:2), plg=list(x=6.2, y=50.2, ncol=2), main="")
 
plot(v, 1, plg=list(x=5.9, y=49.37, horiz=TRUE, cex=1.1), main="", mar=c(5,2,0.5,0.5))

plot(v, density=1:12, angle=seq(18, 360, 20), col=rainbow(12))

plot(v, "NAME_2", col=rainbow(12), border=c("gray", "blue"), lwd=3, type="classes")

plot(v, "AREA", type="interval", breaks=3, mar=c(3.1, 3.1, 2.1, 3.1),
  plg=list(x="topright"), main="")

plot(v, "AREA", type="interval", breaks=c(0,200,250,350), mar=c(2,2,2,2), 
  plg=list(legend=c("<200", "200-250", ">250"), cex=1, 
  bty="o", x=6.4, y=50.125, box.lwd=2, bg="light yellow", title="My Legend"))
}

\keyword{methods}
\keyword{spatial}


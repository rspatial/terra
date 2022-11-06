\name{plot_vector}

\docType{methods}

\concept{plot}

\alias{plot,SpatVector,missing-method}
\alias{plot,SpatVector,numeric-method}
\alias{plot,SpatVector,data.frame-method}
\alias{plot,SpatVector,character-method}
\alias{plot,SpatVectorProxy,missing-method}

\title{Make a map (SpatVector)}

\description{
Plot the values of a SpatVector to make a map. 

See \code{\link{points}}, \code{\link{lines}} or \code{\link{polys}} to add a SpatVector to an existing map.

See \code{\link{plot}} for plotting other object types.
}

\usage{
\S4method{plot}{SpatVector,character}(x, y, col=NULL, type, mar=NULL, legend=TRUE, add=FALSE, axes=!add, 
    main=y, buffer=TRUE, background=NULL, grid=FALSE, ext=NULL, sort=TRUE, 
	decreasing=FALSE, plg=list(), pax=list(), nr, nc, colNA=NA, alpha=NULL, ...)

\S4method{plot}{SpatVector,numeric}(x, y, ...)

\S4method{plot}{SpatVector,missing}(x, y, values=NULL, ...)
}

\arguments{
  \item{x}{SpatVector}
  \item{y}{missing or positive integer or name indicating the layer(s) to be plotted}
  \item{col}{character. Colors. The default for polygon areas is \code{grDevices::rainbow(100, start=.1, end=0.9)}}
  \item{type}{character. Type of map/legend. One of "continuous", "classes", or "interval". If not specified, the type is chosen based on the data}
  \item{mar}{numeric vector of length 4 to set the margins of the plot (to make space for the legend). The default is (3.1, 3.1, 2.1, 7.1) for a single plot with a legend and (3.1, 3.1, 2.1, 2.1) otherwise. Use \code{mar=NA} to not set the margins}
  \item{legend}{logical or character. If not \code{FALSE} a legend is drawn. The character value can be used to indicate where the legend is to be drawn. For example "topright" or "bottomleft". Use \code{plg} for more refined placement. Not supported for continuous legends (the default for raster data)}
  \item{axes}{logical. Draw axes?}
  \item{buffer}{logical. If \code{TRUE} the plotting area is slightly larger than the extent of \code{x}}
  \item{background}{background color. Default is no color (white)}
  \item{plg}{list with parameters for drawing the legend. See the arguments for \code{\link{legend}}. Argument \code{digits} can be used to set the number of digits to print (after the decimal point for a continuous legend)}
  \item{pax}{list with parameters for drawing axes. See the arguments for \code{\link{axis}}. Arguments \code{side}, \code{tick} and \code{lab} can be used to indicate for which of the four axes to draw a line (side), tick mark, or label. The default is \code{c(1:4)} for side and \code{1:2} for the other two. If \code{side} is changed the other two default to that value. Logical argument \code{retro} can be used to use a sexagesimal notation for the labels (degrees/minutes/hemisphere) instead of the standard decimal notation}  
  \item{colNA}{character. color for the NA values}
  \item{alpha}{a single numeric between 0 and 1 to set the transparency for all colors (0 is transparent, 1 is opaque)}
  \item{sort}{logical. If \code{TRUE} legends with categorical values are sorted. You can also supply a vector of the unique values, in the order in which you want them to appear in the legend}
  \item{decreasing}{logical. If \code{TRUE}, legends are sorted in decreasing order}
  \item{grid}{logical. If \code{TRUE} grid lines are drawn. Their properties such as type and color can be set with the \code{pax} argument}
  \item{nc}{positive integer. Optional. The number of columns to divide the plotting device in (when plotting multiple layers)}
  \item{nr}{positive integer. Optional. The number of rows to divide the plotting device in (when plotting multiple layers)}
  \item{main}{character. Main plot titles (one for each layer to be plotted)}
  \item{add}{logical. If \code{TRUE} add the object to the current plot}
  \item{ext}{SpatExtent. Can be use instead of xlim and ylim to set the extent of the plot}
  \item{values}{Either a vector with values to be used for plotting or a two-column data.frame, where the first column matches a variable in \code{x} and the second column has the values to be plotted}
  \item{...}{additional graphical arguments}  
}

\seealso{
\code{\link{points}}, \code{\link{lines}}, \code{\link{polys}}, \code{\link{plot}},

scale bar: \code{\link{sbar}}, north arrow: \code{\link{north}}
}

\examples{ 
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)

plot(v)

plot(v, 2, pax=list(side=1:2), plg=list(x=6.2, y=50.2, cex=1.2))

plot(v, 4, pax=list(side=1:2), plg=list(x=6.2, y=50.2, ncol=2), main="")
 
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


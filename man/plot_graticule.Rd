\name{plot_graticule}

\docType{methods}

\alias{plot,SpatGraticule,missing-method}
\alias{lines,SpatGraticule,missing-method}

\title{Plot a graticule}

\description{
Plot a SpatGraticule. You can create a SpatGraticule with \code{\link{graticule}}.
}

\usage{
\S4method{plot}{SpatGraticule,missing}(x, y, background=NULL, col="black", mar=NULL, labels=TRUE,
	retro=FALSE, lab.loc=c(1,1), lab.lon=NULL, lab.lat=NULL, lab.cex=0.65, 
	lab.col="black", off.lat=0.25, off.lon=0.25, box=FALSE, box.col="black",
	add=FALSE, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{y}{missing or positive integer or name indicating the layer(s) to be plotted}  
  \item{background}{background color. If NULL, no background is drawn}
  \item{mar}{numeric vector of length 4 to set the margins of the plot. To make space for the legend you may use something like \code{c(3.1, 3.1, 2.1, 7.1)}. To fill the plotting canvas, you can use \code{c(0,0,0,0}. Use \code{NA} to not set the margins}
  \item{col}{character. Color for the graticule lines} 
  \item{labels}{logical. If \code{TRUE}, show graticule labels}
  \item{retro}{logical. If \code{TRUE}, show "retro" instead of decimal labels with the graticule}
  \item{lab.loc}{numeric. The first number indicates where the longitude graticule labels should be drawn (1=bottom, 2=top, NA=not drawn, any other number=top and bottom). The second number indicates where the latitude graticule labels should be drawn (1=left, 2=right, NA=not drawn, any other number=left and right)}	
  \item{lab.lon}{positive integers between 1 and the number of labels, indicating which longitude graticule labels should be included}
  \item{lab.lat}{positive integers between 1 and the number of labels, indicating which latitude graticule labels should be included}
  \item{lab.cex}{double. size of the label font}
  \item{lab.col}{character. color of the labels}
  \item{off.lon}{numeric. longitude labels offset}
  \item{off.lat}{numeric. latitude labels offset}
  \item{box}{logical. If \code{TRUE}, the outer lines of the graticule are drawn on top with a sold line \code{lty=1}}
  \item{box.col}{character. color of the outer lines of the graticule if \code{box=TRUE}}
  \item{add}{logical. Add the graticule to the current plot?}
  \item{...}{additional graphical arguments passed to \code{\link{lines}}}
}

\seealso{
 \code{\link{graticule}}, \code{\link{plot}}, 
 \code{\link{points}}, \code{\link{lines}}, \code{\link{polys}}, \code{\link{image}}, \code{scatter\link[terra:scatter]{plot}}, 
 scale bar: \code{\link{sbar}}, north arrow: \code{\link{north}}

}

\examples{
g <- graticule(60, 30, crs="+proj=robin")

plot(g, background="azure", col="red", lty=2, box=TRUE)
plot(g, background="azure", col="light gray", lab.loc=c(1,2), 
		lab.lon=c(2,4,6), lab.lat=3:5, lty=3, retro=TRUE)
}


\keyword{methods}
\keyword{spatial}


\name{plot}

\docType{methods}

\alias{plot}
\alias{plot,SpatRaster,numeric-method}
\alias{plot,SpatRaster,missing-method}

\title{Plot a SpatRaster}

\description{
Plot the values of a SpatRaster to make a map
}

\usage{
\S4method{plot}{SpatRaster,numeric}(x, y=1, col, type, mar=NULL, legend=TRUE, axes=TRUE, plg=list(), 
     pax=list(), maxcell=50000, smooth=FALSE, range=NULL, levels=NULL, fun=NULL, 
     colNA=NULL, alpha=NULL, ...)

\S4method{plot}{SpatRaster,missing}(x, y, maxcell=50000, main, mar=NULL, nc, nr, maxnl=16, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{missing or positive integer or name indicating the layer(s) to be plotted}
  \item{col}{character. Colors}
  \item{type}{character. Type of map/legend. One of "continuous", "classes", or "interval"}
  \item{mar}{numeric vector of lenght 4 to set the margins of the plot (to make space for the legend). The default is (3.1, 3.1, 2.1, 7.1) for a single plot with a legend and (3.1, 3.1, 2.1, 2.1) otherwise. When multiple layers are plotted, the defaults are 
  Use \code{mar=NA} to not set the margins}
  \item{legend}{logical. Draw a legend?}
  \item{axes}{logical. Draw axes?}
  \item{plg}{list with parameters for drawing the legend. See the arguments for \code{\link{legend}}}
  \item{pax}{list with parameters for drawing axes. See the arguments for \code{\link{axis}}}  
  \item{maxcell}{positive integer. Maximum number of cells to use for the plot}
  \item{smooth}{logical. If \code{TRUE} the cell values are smoothed (for continuous legend)}
  \item{range}{numeric. minimum and maximum values to be used for the continuous legend }
  \item{levels}{character. labels to be used for the classes legend }
  \item{fun}{function to be called after plotting each SpatRaster layer to add something to each map (such as text, legend, lines). For example, with SpatVector \code{v}, you could do \code{fun=function() lines(v)}. The function may have one argument, representing the the layer that is plotted (1 to the number of layers) }  
  \item{colNA}{character. color for the NA values}
  \item{alpha}{numeric between 0 and 1 to set the transparency for all colors (0 is transparent, 1 is opaque)}
  \item{nc}{positive integer. Optional. The number of columns to divide the plotting device in (when plotting multiple layers)}
  \item{nr}{positive integer. Optional. The number of rows to divide the plotting device in (when plotting multiple layers)}
  \item{main}{character. Main plot titles (one for each layer to be plotted)}
  \item{maxnl}{positive integer. Maximum number of layers to plot (for a multi-layer object)}
  \item{...}{arguments passed to \code{plot("SpatRaster", "numeric")} and additional graphical arguments}  
}

\seealso{
 \code{\link{image}}, \code{\link[terra:plotVector]{plotVector}}, \code{scatter\link[terra:scatter]{plot}}, \code{\link{sbar}}
}

\examples{ 
f <- system.file("ex/elev.tif", package="terra") 
r <- rast(f)
plot(r)

plot(r, type="interval")

e <- c(6.3, 6.35, 49.9, 50.1)
plot(r, plg=list(ext=e, title="Title\n", title.cex=1.25), pax=list(sides=1:2))


d <- classify(r, c(100,200,300,400,500,600))
plot(d, type="classes")

plot(d, type="interval", breaks=1:5) 
plot(d, type="interval", breaks=c(1,4,5), plg=list(legend=c("1-4", "4-5"))) 
plot(d, type="classes", plg=list(legend=c("Mr", "Xx", "As", "Zx", "Bb"), x="bottomright"))

x <- trunc(r/200)
x <- as.factor(x)
levels(x) <- c("earth", "wind", "fire")
plot(x, plg=list(x="topright"),mar=c(2,2,2,2))

# two plots with the same legend
dev.new(width=6, height=4, noRStudioGD = TRUE)
par(mfrow=c(1,2))
plot(r, range=c(50,600))
plot(r/2, range=c(50,600))

# as you only need one legend:
par(mfrow=c(1,2))
plot(r, range=c(50,600), mar=c(4, 3, 4, 3), plg=list(shrink=0.9, cex=.8), 
	pax=list(sides=1:2, cex.axis=.6))
#text(182500, 335000, "Two maps, one plot", xpd=NA)
plot(r/2, range=c(50,600), mar=c(4, 2, 4, 4), legend=FALSE, 
	pax=list(sides=c(1,4), cex.axis=.6))
}

\keyword{methods}
\keyword{spatial}


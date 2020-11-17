\name{plot}

\docType{methods}

\alias{plot}
\alias{plot,SpatRaster,missing-method}
\alias{plot,SpatRaster,numeric-method}


\title{Plot a SpatRaster}

\description{
Plot the values of a SpatRaster to make a map
}

\usage{
\S4method{plot}{SpatRaster,numeric}(x, y=1, col, type, mar=c(5.1, 4.1, 4.1, 7.1), legend=TRUE, axes=TRUE,
  pal=list(), pax=list(), maxcell=50000, smooth=FALSE, range=NULL, levels=NULL, fun=NULL, ...)

\S4method{plot}{SpatRaster,missing}(x, y, col, type, mar=c(5.1, 4.1, 4.1, 7.1), legend=TRUE, axes=TRUE,
  pal=list(), pax=list(), maxcell=50000, smooth=FALSE, range=NULL, levels=NULL, fun=NULL, nc, nr, main, maxnl=16, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{missing or positive integer indicating the layer(s) to be plotted, or the name of the layer to be mapped}
  \item{col}{character. Colors}
  \item{type}{character. Type of map/legend. One of "continuous", "classes", or "interval"}
  \item{mar}{numeric vector of lenght 4 to set the margins of the plot}
  \item{legend}{logical. Draw a legend?}
  \item{axes}{logical. Draw axes?}
  \item{pal}{list. Parameters for the legend}
  \item{pax}{list. Parameters for drawing axes}
  \item{maxcell}{positive integer. Maximum number of cells to use for the plot}
  \item{smooth}{logical. If \code{TRUE} the cell values are smoothed (for continuous legend)}
  \item{range}{numeric. minimum and maximum values to be used for the continuous legend }
  \item{levels}{character. labels to be used for the classes legend }
  \item{fun}{function to be called after plotting each SpatRaster layer to add something to each map (such as text, legend, lines). For example, with SpatVector \code{v}, you could do \code{fun=function() lines(v)}. The function may have one argument, representing the the layer that is plotted (1 to the number of layers) }  
  \item{nc}{positive integer. Optional. The number of columns to divide the plotting device in (when plotting multiple layers)}
  \item{nr}{positive integer. Optional. The number of rows to divide the plotting device in (when plotting multiple layers)}
  \item{main}{character. Main plot titles (one for each layer to be plotted)}
  \item{maxnl}{positive integer. Maximum number of layers to plot (for a multi-layer object)}
  \item{...}{additional graphical arguments}  
}

\seealso{
 \code{\link{image}}, \code{\link[terra:plotVector]{plotVector}}, \code{scatter\link[terra:scatter]{plot}}
}



\examples{ 
f <- system.file("ex/test.tif", package="terra") 
r <- rast(f)
plot(r)

plot(r, type="interval")

e <- c(178200,178400,331600,333600)
plot(r, pal=list(ext=e, title="Title\n", title.cex=1.5))

d <- (r > 400) + (r > 600)
plot(d)
plot(d, type="classes")

plot(d, type="interval", breaks=0:3) 
plot(d, type="interval", breaks=3, pal=list(legend=c("0-1", "1-2", "2-3"))) 
plot(d, type="classes", pal=list(legend=c("M", "X", "A")))

x <- trunc(r/600)
x <- as.factor(x)
levels(x) <- c("earth", "wind", "fire")
plot(x)


# two plots with the same legend
dev.new(width=6, height=4, noRStudioGD = TRUE)
par(mfrow=c(1,2))
plot(r, range=c(100,1800))
plot(r/2, range=c(100,1800))

# as you only need one legend:
par(mfrow=c(1,2))
plot(r, range=c(100,1800), mar=c(4, 3, 4, 3), pal=list(shrink=0.9, cex=.8), 
	pax=list(sides=1:2, cex.axis=.6))
text(182500, 335000, "Two maps, one plot", xpd=NA)
plot(r/2, range=c(100,1800), mar=c(4, 2, 4, 4), legend=FALSE, 
	pax=list(sides=c(1,4), cex.axis=.6))
}

\keyword{methods}
\keyword{spatial}


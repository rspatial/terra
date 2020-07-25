\name{plotRaster}

\docType{methods}

\alias{plot}
\alias{plot,SpatRaster,missing-method}
\alias{plot,SpatRaster,numeric-method}

\alias{image}
\alias{image,SpatRaster-method}


\title{Plot a SpatRaster}

\description{
Plot (that is, make a map of) the values of a SpatRaster or SpatVector
}

\usage{
\S4method{plot}{SpatRaster,numeric}(x, y=1, col, 
  type, mar=c(5.1, 4.1, 4.1, 7.1), legend=TRUE, axes=TRUE,
  pal=list(), pax=list(), maxcell=50000, ...)

\S4method{plot}{SpatRaster,missing}(x, y, maxcell=50000, nc, nr, main, 
   maxnl=16, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{missing or positive integer indicating the layer(s) to be plotted, or the name of the layer to be mapped}
  \item{col}{character. Colors}
  \item{type}{character. Type of map/legend. One of "continuous", "classes", or "internval"}
  \item{mar}{numeric vector of lenght 4 to set the margins of the plot}

  \item{legend}{logical. Draw a legend?}
  \item{axes}{logical. Draw axes?}
  \item{pal}{list. Parameters for the legend}
  \item{pax}{list. Parameters for drawing axes}
  
   \item{maxcell}{positive integer. Maximum number of cells to use for the plot}
   
  \item{nc}{positive integer. Optional. The number of columns to divide the plotting device in (when plotting multiple layers)}
  \item{nr}{positive integer. Optional. The number of rows to divide the plotting device in (when plotting multiple layers)}
  \item{main}{character. Main plot titles (one for each layer to be plotted)}
  \item{maxnl}{positive integer. Maximum number of layers to plot (for a multi-layer object)}


  \item{...}{additional graphical arguments}  
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


plot(d, type="interval", levels=0:3) 
plot(d, type="interval", levels=3, pal=list(legend=c("0-1", "1-2", "2-3"))) 
plot(d, type="classes", pal=list(legend=c("M", "X", "A")))

x <- trunc(r/600)
x <- as.factor(x)
levels(x) <- c("earth", "wind", "fire")
plot(x)

}

\keyword{methods}
\keyword{spatial}


\name{as.raster}

\alias{as.raster}
\alias{as.raster,SpatRaster-method}
  
\title{Coerce to a "raster" object}

\description{
Implementation of the generic \code{\link[grDevices]{as.raster}} function to create a "raster" (small r) object. Such objects can be used for plotting with the \code{\link[graphics]{rasterImage}} function. NOT TO BE CONFUSED with the Raster* (big R) objects defined by the 'raster' package! 
}

\usage{
\S4method{as.raster}{SpatRaster}(x, maxcell=500000, col)
}


\arguments{
 \item{x}{ SpatRaster }
 \item{maxcell}{positive integer. Maximum number of cells to use for the plot}
 \item{col}{vector of colors. The default is \code{map.pal("viridis", 100)}}
}

 

\value{
'raster' object
}



\examples{
r <- rast(ncols=3, nrows=3)
values(r) <- 1:ncell(r)
as.raster(r)
}

\keyword{spatial}
\keyword{methods}

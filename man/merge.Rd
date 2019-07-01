\name{merge}

\docType{methods}

\alias{merge}
\alias{merge,SpatRaster,SpatRaster-method}
\alias{merge,SpatExtent,SpatExtent-method}

\title{ 
Merge SpatRaster or SpatExtent objects
}

\description{
Merge SpatRasters to form a new SpatRaster object with a larger spatial extent. If objects overlap, the values get priority in the same order as the arguments. See \code{\link[raster]{reclassify}} to merge a \code{SpatRaster} and a \code{data.frame}. You can also merge SpatExtent objects.
}

\usage{
\S4method{merge}{SpatRaster,SpatRaster}(x, y, ..., filename="", overwrite=FALSE, wopt=list())
\S4method{merge}{SpatExtent,SpatExtent}(x, y, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatExtent object}
  \item{y}{object of same class as \code{x}}
  \item{...}{additional objects of the same class as \code{x}}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
} 

\details{
The SpatRaster objects must have the same origin and spatial resolution. In areas where the SpatRaster objects overlap, the values of the Raster object that is last in the sequence of arguments will be retained. 
}

\value{
SpatRaster or SpatExtent
}


\examples{
x <- rast(xmax=-150, ymin=60, ncols=30, nrows=30)
y <- rast(xmin=-100, xmax=-50, ymax=50, ymin=30)
y <- setRes(y, res(x))
values(x) <- 1:ncell(x)
values(y) <- 1:ncell(y)
mr <- merge(x, y)
plot(mr)

# if you have many SpatRaster objects in a list
# you can use do.call:
s <- list(x, y)
# add arguments such as filename
# s$filename <- 'test.tif'
m <- do.call(merge, s)
}

\keyword{methods}
\keyword{spatial}


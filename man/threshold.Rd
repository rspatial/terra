\name{threshold}

\docType{methods}

\alias{threshold}
\alias{threshold,SpatRaster-method}


\title{Thresholding}

\description{
Compute a threshold to divide the values of a SpatRaster into two groups, and use that threshold to classify the raster. 
}

\usage{
\S4method{threshold}{SpatRaster}(x, method="otsu", maxcell=1000000, combine=FALSE, as.raster=TRUE, filename="", ...)
}


\arguments{
\item{x}{SpatRaster}
\item{method}{character. One of "mean", "median" or "otsu"}
\item{maxcell}{positive integer. Maximum number of cells to use to compute the threshold}
\item{combine}{logical. If \code{TRUE} the layers of \code{x} are combine to compute a single threshold}
\item{as.raster}{logical. If \code{TRUE} a classified SpatRaster is returned. Otherwise the threshold(s) are retured}
\item{filename}{character. Output filename}
\item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}  
}

\value{
numeric or SpatRaster
}


\examples{

s <- rast(system.file("ex/logo.tif", package="terra"))
threshold(s, "mean", as.raster=FALSE)
threshold(s, "mean", combine=TRUE, as.raster=FALSE)


plot(threshold(s, "otsu"))

}

\keyword{spatial}


\name{twoClasses}

\docType{methods}

\alias{twoClasses}
\alias{twoClasses,SpatRaster-method}


\title{Thresholding}

\description{
Compute a threshold to divide the values of a SpatRaster into two groups, and use that threshold to classify the raster. 
}

\usage{
\S4method{twoClasses}{SpatRaster}(x, method="otsu", maxcell=1000000, combine=FALSE, 
			as.raster=TRUE, filename="", ...)
}


\arguments{
\item{x}{SpatRaster}
\item{method}{character. One of "mean", "median" or "otsu" for Otsu's method}
\item{maxcell}{positive integer. Maximum number of cells to use to compute the threshold}
\item{combine}{logical. If \code{TRUE} the layers of \code{x} are combine to compute a single threshold}
\item{as.raster}{logical. If \code{TRUE} a classified SpatRaster is returned. Otherwise the threshold(s) are retured}
\item{filename}{character. Output filename}
\item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}  
}

\value{
numeric or SpatRaster
}

\seealso{\code{\link{divide}}}

\examples{

s <- rast(system.file("ex/logo.tif", package="terra"))
twoClasses(s, "mean", as.raster=FALSE)
twoClasses(s, "mean", combine=TRUE, as.raster=FALSE)


plot(twoClasses(s, "otsu"))

}

\keyword{spatial}


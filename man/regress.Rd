\name{regress}

\docType{methods}

\alias{regress}
\alias{regress,SpatRaster,SpatRaster-method}
\alias{regress,SpatRaster,numeric-method}
\alias{regress,SpatRaster,missing-method}

\title{Cell level regression}

\description{
Run a regression model for each cell of a SpatRaster. The independent variable can either be defined by a vector, or another SpatRaster to make it spatially variable.
} 

\usage{

\S4method{regress}{SpatRaster,numeric}(y, x=1:nlyr(y), formula=y~x, na.rm=FALSE, cores=1, filename="", overwrite=FALSE, ...)

\S4method{regress}{SpatRaster,SpatRaster}(y, x, formula=y~x, na.rm=FALSE, cores=1, filename="", overwrite=FALSE)  
}


\arguments{
\item{y}{SpatRaster}
\item{x}{SpatRaster or numeric. If missing, \code{as.numeric(time(y)} is used or \code{1:nlyr(y)} if \code{y} does not have timestamps}
\item{formula}{regression formula in the general form of \code{y ~ x}. You can add additional terms such as \code{I(x^2)} }
\item{na.rm}{logical. Remove NA values?}
\item{cores}{positive integer. If \code{cores > 1}, a 'parallel' package cluster with that many cores is created and used. You can also supply a cluster object.}
\item{filename}{character. Output filename}
\item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
\item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
}

\value{SpatRaster}


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
x <- regress(s)
}

\keyword{methods}
\keyword{spatial}


\name{sds}

\docType{methods}

\alias{sds}
\alias{sds,character-method}
\alias{sds,list-method}
\alias{sds,SpatRaster-method}
\alias{[<-,SpatDataSet,numeric,missing-method}
\alias{sds_info}


\title{Create a SpatDataSet}

\description{
Methods to create a SpatDataSet. This is an object to hold "sub-datasets", each a SpatRaster that in most cases will have multiple layers. 

sds_info provides information about sub-datasets in a file.

}

\usage{
\S4method{sds}{character}(x, ids=0, ...)

\S4method{sds}{SpatRaster}(x, ...) 

\S4method{sds}{list}(x, ...) 

sds_info(filename, ...)
}

\arguments{
\item{x}{character (filename) or SpatRaster, or list of SpatRaster objects}
\item{ids}{optional. vector of integer subdataset ids. Ignored if the first value is not a postive integer}
\item{...}{additional arguments. Can be other \code{SpatRaster} objects if \code{x} is a \code{SpatRaster}} 
\item{filename}{character}
}


\value{
SpatDataSet
}


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
x <- sds(s, s/2)
names(x) <- c("first", "second")
x
length(x)

# extract the second SpatRaster
x[2]
}


\keyword{methods}
\keyword{spatial}


\name{rstk}

\docType{methods}

\alias{rstk}
\alias{rstk,character-method}
\alias{rstk,list-method}
\alias{rstk,SpatRaster-method}
\alias{[<-,SpatStack,numeric,missing-method}
\alias{sds_info}


\title{Create a SpatStack}

\description{
Methods to create a SpatStack. This is an object to hold "sub-datasets", each a SpatRaster that in most cases will have multiple layers. 

sds_info provides information about sub-datasets in a file.

}

\usage{
\S4method{rstk}{character}(x, ids=0, ...)

\S4method{rstk}{SpatRaster}(x, ...) 

\S4method{rstk}{list}(x, ...) 

sds_info(filename, ...)
}

\arguments{
\item{x}{character (filename) or SpatRaster, or list of SpatRaster objects}
\item{ids}{optional. vector of integer subdataset ids. Ignored if the first value is not a postive integer}
\item{...}{additional arguments. Can be other \code{SpatRaster} objects if \code{x} is a \code{SpatRaster}} 
\item{filename}{character}
}


\value{
SpatStack
}


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
x <- rstk(s)
names(x) <- "first"
y <- c(x, second=s)

# extract the second SpatRaster
y[2]
}


\keyword{methods}
\keyword{spatial}


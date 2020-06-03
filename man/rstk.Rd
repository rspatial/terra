\name{rstk}

\docType{methods}

\alias{rstk}
\alias{rstk,character-method}
\alias{rstk,SpatRaster-method}
\alias{[<-,SpatStack,numeric,missing-method}

\title{Create a SpatStack}

\description{
Methods to create a SpatStack. This is an object to hold "sub-datasets", each a SpatRaster that in most cases will have multiple layers. 
}

\usage{
\S4method{rstk}{character}(x, ...)

\S4method{rstk}{SpatRaster}(x, name="sd1", ...) 
}

\arguments{
\item{x}{character (filename) or SpatRaster}
\item{name}{character to name the sub-dataset}
\item{...}{additional arguments. Can be other \code{SpatRaster} objects if \code{x} is a \code{SpatRaster}}  
}


\value{
SpatStack
}


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
x <- rstk(s, "first")
y <- c(x, second=s)

# extract the second SpatRaster
y[2]
}


\keyword{methods}
\keyword{spatial}


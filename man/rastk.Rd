\name{rast}

\docType{methods}

\alias{rastk}
\alias{rastk,character-method}
\alias{rastk,SpatRaster-method}

\title{Create a SpatRasterStack}

\description{
Methods to create a SpatRasterStack. This is an object to hold "sub-datasets", each a SpatRaster that in most cases will have multiple layers. 
}

\usage{
\S4method{rastk}{character}(x, ...)

\S4method{rastk}{SpatRaster}(x, ...) 
}

\arguments{
\item{x}{filename (character) or SpatRaster}
\item{...}{additional arguments. None implemented }  
}


\value{
SpatRasterStack
}


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
x <- rastk(s, "first")
add(x, s, "second")
}


\keyword{methods}
\keyword{spatial}


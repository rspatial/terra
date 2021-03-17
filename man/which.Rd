\name{which}

\docType{methods}

\alias{Which}
\alias{Which,SpatRaster-method}

\title{Which cells are TRUE?}

\description{
\code{Which} returns a single layer SpatRaster with cell values that are the first layer in the input that is not 0 (are \code{TRUE}) 
}


\usage{
\S4method{Which}{SpatRasterLayer}(x)
}

\arguments{
  \item{x}{SpatRaster}
}


\value{
SpatRaster
}



\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
x <- Which(s)
}

\keyword{spatial}

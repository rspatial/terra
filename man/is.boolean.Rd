\name{is.boolean}

\docType{methods}

\alias{is.boolean}
\alias{is.boolean,SpatRaster-method}
\alias{as.boolean}
\alias{as.boolean,SpatRaster-method}

  
\title{Logical rasters}

\description{

A SpatRaster layer can reprsent a logical (Boolean) variable}

\usage{
\S4method{is.boolean}{SpatRaster}(x)
}


\arguments{
  \item{x}{SpatRaster}
}  


\value{
logical
}

\examples{
r <- rast(nrows=10, ncols=10, vals=1:100)
is.boolean(r)

x <- r > 25
is.boolean(x)
}

\keyword{methods}
\keyword{spatial}


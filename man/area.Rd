\name{diff}

\alias{diff}
\alias{diff,SpatRaster-method}
  
\title{diff}

\description{
Returns lagged differences.
}

\usage{
\S4method{diff}{SpatRaster}(x, filename="", overwrite=FALSE,  wopt=list(), ...)
}


\arguments{
  \item{x}{SpatRaster}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented} 
}

\value{
SpatRaster
}



\examples{
s <- rast(system.file("exdata/logo.tif", package="terra"))   
d <- diff(s)
}

\keyword{methods}
\keyword{spatial}

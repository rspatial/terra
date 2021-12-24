\name{fillTime}

\docType{methods}

\alias{fillTime}
\alias{fillTime,SpatRaster-method}

\title{ 
Fill time gaps in a SpatRaster
}

\description{
Add empty layers in between existing layers such that the time step between each layer is the same. 
See \code{\link{approximate}} to estimate values for these layer (and other missing values)
}

\usage{
\S4method{fillTime}{SpatRaster}(x, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{filename}{character. Output filename}
  \item{...}{list with named options for writing files as in \code{\link{writeRaster}}}
} 

\value{
SpatRaster
}

\seealso{
\code{\link{approximate}}
}

\examples{
r <- rast(system.file("ex/logo.tif", package="terra"))   
s <- c(r, r)
time(s) <- as.Date("2001-01-01") + c(0:2, 5:7)
time(s)
ss <- fillTime(s)
time(ss)

a <- approximate(ss)
}

\keyword{methods}
\keyword{spatial}


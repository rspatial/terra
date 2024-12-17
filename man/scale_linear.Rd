\name{scale_linear}

\alias{scale_linear}
\alias{scale_linear,SpatRaster-method}


\title{Scale values linearly}

\description{
Linear scaling of raster cell values between a specified minimum and maximum value.
}

\usage{
\S4method{scale_linear}{SpatRaster}(x, min=0, max=1, filename="", ...)
}


\arguments{
 \item{x}{SpatRaster}
 \item{min}{minimum value to scale to}
 \item{max}{maximum value to scale to}
 \item{filename}{character. Output filename}
 \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}
 
\value{
SpatRaster
}


\seealso{ \code{\link{scale}} }

\examples{
r <- rast(system.file("ex/logo.tif", package="terra"))   
s1 <- scale_linear(r)
s2 <- scale_linear(r, 1, 10)
}

\keyword{ spatial }

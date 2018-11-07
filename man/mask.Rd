\name{mask}

\docType{methods}

\alias{mask}
\alias{mask,SpatRaster,SpatRaster-method}

\title{Mask values in a SpatRaster object}

\description{
Create a new SpatRaster that has the same values as SpatRaster \code{x}, except for the cells that are \code{NA} in another SpatRaster (the 'mask'). These cells become \code{NA}.
} 

\usage{
\S4method{mask}{SpatRaster,SpatRaster}(x, mask, filename="", overwrite=FALSE, ...)
}

\arguments{
\item{x}{SpatRaster object}
\item{mask}{SpatRaster object}
\item{filename}{character. Optional output filename}
\item{overwrite}{Logical. If \code{TRUE}, \code{filename} is overwritten}
\item{...}{additional arguments as in \code{\link{writeRaster}}}
}

\value{SpatRaster object}

\seealso{\link{crop}}


\examples{
r <- rast(ncol=10, nrow=10)
m <- rast(ncol=10, nrow=10)
values(r) <- 1:100
x <- runif(ncell(r))
x[x < 0.5] <- NA
values(m) <- x
mr <- mask(r, m)
}

\keyword{methods}
\keyword{spatial}


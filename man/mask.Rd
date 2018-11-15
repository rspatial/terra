\name{mask}

\docType{methods}

\alias{mask}
\alias{mask,SpatRaster,SpatRaster-method}

\title{Mask values in a SpatRaster}

\description{
Create a new SpatRaster that has the same values as SpatRaster \code{x}, except for the cells that are \code{NA} in another SpatRaster (the 'mask'). These cells become \code{NA}.
} 

\usage{
\S4method{mask}{SpatRaster,SpatRaster}(x, mask, filename="", overwrite=FALSE, writeopt, ...)
}

\arguments{
\item{x}{SpatRaster}
\item{mask}{SpatRaster}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{writeopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{Additional arguments. None implemented}  

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


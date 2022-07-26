\name{concats}

\docType{methods}

\alias{concats}
\alias{concats,SpatRaster-method}
  
  
\title{Concatenate categorical rasters}

\description{
Combine two categorical rasters by concatenating their levels.
}

\usage{
\S4method{concats}{SpatRaster}(x, y, filename="", ...)
}


\arguments{
  \item{x}{SpatRaster (with a single, categorical, layer)}
  \item{y}{SpatRaster (with a single, categorical, layer)}
  \item{filename}{character. Output filename}  
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}  
}  


\value{
SpatRaster
}

\seealso{\code{\link{cats}}}


\examples{
set.seed(0)
r <- rast(nrows=10, ncols=10)
values(r) <- sample(3, ncell(r), replace=TRUE)
levels(r) <- data.frame(id=1:3, cover=c("forest", "water", "urban"))

rr <- rast(r)
values(rr) <- sample(1:3, ncell(rr), replace=TRUE)
levels(rr) <- data.frame(id=c(1:3), color=c("red", "green", "blue"))

x <- concats(r, rr)
x
levels(x)[[1]]
}

\keyword{methods}
\keyword{spatial}


\name{xapp}

\docType{methods}

\alias{xapp}
\alias{xapp,SpatRaster,SpatRaster-method}

\title{Apply a function to the cells of a two SpatRasters}

\description{ 
Apply a function to the values of each cell of two (multilayer) SpatRasters. 
}

\usage{
\S4method{xapp}{SpatRaster,SpatRaster}(x, y, fun, ..., filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatRaster with the same geometry as \code{x}}
  \item{fun}{a function that operates on two vectors}
  \item{...}{additional arguments for \code{fun}. These are typically numerical constants. They should *never* be another SpatRaster}  
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}


\seealso{\code{\link{app}}, \code{\link{lapp}}, \code{\link{tapp}}, \code{\link[terra]{Math-methods}}, \code{\link{roll}} }


\examples{
r <- rast(ncols=10, nrows=10, nlyr=5)
set.seed(1)
r <- init(r, runif)
s <- init(r, runif)
x <- xapp(r, s, fun=cor)
}


\keyword{methods}
\keyword{spatial}

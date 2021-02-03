\name{arith}

\docType{methods}

\alias{arith}
\alias{arith,SpatRaster-method}

\title{apply a local function}

\description{ 
Apply a function that (arithmetically) operates to individual cells and layers of a SpatRaster, to return the same number of layers as in the input SpatRaster. That, is unlike with \code{\link{app}}, without summarizing over layers.
}

\usage{
\S4method{arith}{SpatRaster}(x, fun, ..., filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{fun}{function}
  \item{...}{additional arguments for \code{fun}}  
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{ \code{\link{arith}} }


\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)

x <- arith(r, sqrt)

# equivalent to 
y <- sqrt(x)
}


\keyword{methods}
\keyword{spatial}

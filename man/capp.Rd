\name{capp}

\docType{methods}

\alias{capp}
\alias{capp,SpatRaster-method}

\title{apply a local function}

\description{ 
Apply a function operates to individual cells and layers of a SpatRaster, to return the same number of layers as in the input SpatRaster. That, is unlike with \code{\link{app}}, without summarizing over layers.
}

\usage{
\S4method{capp}{SpatRaster}(x, fun, ..., filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{fun}{function}
  \item{...}{additional arguments for \code{fun}}  
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{ \code{\link{app}} }


\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)

x <- capp(r, sqrt)

# equivalent to 
y <- sqrt(x)
}


\keyword{methods}
\keyword{spatial}

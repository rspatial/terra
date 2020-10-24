\name{loca}

\docType{methods}

\alias{loca}
\alias{loca,SpatRaster-method}

\title{apply a local function}

\description{ 
A "local" function operates on individual cells and layers of a SpatRaster, thus the number of layers returned is the same as the input (unlike for \code{\link{app}})

The method is not called \code{local} to avoiding hiding the base function with that name.
}

\usage{
\S4method{loca}{SpatRaster}(x, fun, ..., filename="", overwrite=FALSE, wopt=list())
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

x <- loca(r, sqrt)

# equivalent to 
y <- sqrt(x)
}


\keyword{methods}
\keyword{spatial}

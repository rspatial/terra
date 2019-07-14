\name{app}

\docType{methods}

\alias{app}
\alias{app,SpatRaster-method}

\title{Apply a function to the cells of a SpatRaster}

\description{ 
Apply a function to values of a SpatRaster. Similar to \code{\link[base]{apply}} 
}

\usage{
\S4method{app}{SpatRaster}(x, fun, ..., filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster object}
  \item{fun}{function}
  \item{...}{additional arguments for \code{fun}}  
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{ \code{ \link[terra]{overlay}}, \link[terra]{math} }


\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
x <- c(r, r, r)
s <- app(x, fun=sum)
}


\keyword{methods}
\keyword{spatial}

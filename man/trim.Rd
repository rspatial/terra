\name{trim}

\alias{trim}
\alias{trim,SpatRaster-method}

\title{Trim a SpatRaster}

\description{
Trim (shrink) a SpatRaster by removing outer rows and columns that are \code{NA} or another value.
}

\usage{
\S4method{trim}{SpatRaster}(x, padding=0, value=NA, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{padding}{integer. Number of outer rows/columns to keep}
  \item{value}{numeric. The value of outer rows or columns that are to be removed}
   \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}


\value{
SpatRaster
}


\examples{
r <- rast(ncols=10, nrows=10, xmin=0,xmax=10,ymin=0,ymax=10)
v <- rep(NA, ncell(r))
v[c(12,34,69)] <- 1:3
values(r) <- v
s <- trim(r) 
}

\keyword{spatial}


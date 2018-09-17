\name{trim}

\alias{trim}
\alias{trim,SpatRaster-method}

\title{Trim a SpatRaster}

\description{
Trim (shrink) a SpatRaster object by removing outer rows and columns that are \code{NA}.
}

\usage{
\S4method{trim}{SpatRaster}(x, padding=0, filename='', ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{padding}{integer. Number of outer rows/columns to keep}
  \item{filename}{character. Optional output filename}
  \item{...}{additional arguments as for \code{\link{writeRaster}}}
}


\value{
SpatRaster
}


\examples{
r <- rast(ncol=18,nrow=18)
v <- rep(NA, ncell(r))
v[39:40] <- 1
v[205] <- 2
values(r) <- v
s <- trim(r) 
}

\keyword{spatial}


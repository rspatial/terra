\name{trim}

\alias{trim}
\alias{trim,SpatRaster-method}

\title{Trim a SpatRaster}

\description{
Trim (shrink) a SpatRaster object by removing outer rows and columns that are \code{NA}.
}

\usage{
\S4method{trim}{SpatRaster}(x, padding=0, filename="", format="", datatype="FLT4S", overwrite=FALSE, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{padding}{integer. Number of outer rows/columns to keep}
  \item{filename}{Character. Output filename. Optional}
  \item{format}{Character. Output file format. See \code{\link{writeRaster}}}
  \item{datatype}{Character. Output file data type. See \code{\link{writeRaster}}} 
  \item{overwrite}{Logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{...}{Additional arguments. None implemented}
}


\value{
SpatRaster
}


\examples{
r <- rast(ncol=10, nrow=10, xmin=0,xmax=10,ymin=0,ymax=10)
v <- rep(NA, ncell(r))
v[c(12,34,69)] <- 1:3
values(r) <- v
s <- trim(r) 
}

\keyword{spatial}


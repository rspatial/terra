\name{focalValues}

\docType{methods}

\alias{focalValues}
\alias{focalValues,SpatRaster-method}
  
\title{Get focal values}

\description{
Get a matrix in which each row had the focal values of a cell. These are the values of a cell and a rectangular window around it.
}

\usage{
\S4method{focalValues}{SpatRaster}(x, w=3, row=1, nrows=nrow(x), fill=NA)
}

\arguments{
 \item{x}{SpatRaster or SpatVector}
 \item{w}{window. The window can be defined as one (for a square) or two odd numbers (row, col); or with an odd sized matrix}
 \item{row}{positive integer. Row number to start from, should be between 1 and nrow(x)}
 \item{nrows}{positive integer. How many rows?}
 \item{fill}{numeric used as values for imaginary cells outside the raster}
}

\value{
matrix
}


\examples{
r <- rast(ncol=4, nrow=4, crs="+proj=utm +zone=1 +datum=WGS84")
values(r) <- 1:ncell(r)
focalValues(r)
}

\keyword{spatial}
\keyword{methods}

\name{cells}

\docType{methods}

\alias{cells}
\alias{cells,SpatRaster,missing-method}
\alias{cells,SpatRaster,SpatVector-method}
\alias{cells,SpatRaster,SpatExtent-method}
\alias{cells,SpatRaster,matrix-method}


\title{Get cell numbers}

\description{ 
Get the cell numbers covered by a SpatVector or SpatExtent. Or all non \code{NA} values.
}

\usage{
\S4method{cells}{SpatRaster,missing}(x, y, ...)

\S4method{cells}{SpatRaster,SpatVector}(x, y, ...)

\S4method{cells}{SpatRaster,matrix}(x, y, ...)

\S4method{cells}{SpatRaster,SpatExtent}(x, y, ...)
}


\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatVector, SpatExtent, 2-column matrix representing points, or missing}
  \item{...}{additional arguments. None implemented}  
}

\value{
matrix
}


\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
r[c(1:25, 31:100)] <- NA
cells(r)

m <- cbind(x=c(0,10,-30), y=c(40,-10,20))
cellFromXY(r, m)

v <- vect(m)
cells(r, v)
}


\keyword{methods}
\keyword{spatial}

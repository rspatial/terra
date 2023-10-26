\name{Arith-methods}

\docType{methods}

\alias{Arith-methods}

\alias{Arith,SpatExtent,numeric-method}
\alias{Arith,numeric,SpatExtent-method}
\alias{Arith,SpatExtent,SpatExtent-method}

\alias{Arith,SpatRaster,SpatRaster-method}
\alias{Arith,SpatVector,SpatVector-method}
\alias{Arith,SpatRaster,missing-method}
\alias{Arith,numeric,SpatRaster-method}
\alias{Arith,SpatRaster,numeric-method}
\alias{Arith,logical,SpatRaster-method}
\alias{Arith,SpatRaster,logical-method}
\alias{Arith,missing,SpatRaster-method}
\alias{Arith,SpatRaster,matrix-method}
\alias{Arith,matrix,SpatRaster-method}


\title{Arithmetic}

\description{

Standard arithmetic operators for computations with SpatRasters. Computations are local (applied on a cell by cell basis). If multiple SpatRaster objects are used, these must have the same geometry (extent and resolution). These operators have been implemented:

\code{ +, -, *, /, ^, \%\%, \%/\% }

You can also use a SpatRaster and a vector or a matrix. If you use a SpatRaster with a vector of multiple numbers, each element in the vector is considered a layer (with a constant value). If you use a SpatRaster with a matrix, the number of columns of the matrix must match the number of layers of the SpatRaster. The rows are used to match the cells. That is, if there are two rows, these match cells 1 and 2, and they are recycled to 3 and 4, etc. 

The following methods have been implemented for \code{(SpatExtent, SpatExtent)}: \code{+, -}, and the following for \code{(SpatExtent, numeric)}: \code{+, -, *, /, \%\%}
}

\seealso{
\code{\link{ifel}} to conveniently combine operations and \code{\link{Math-methods}} or \code{\link{app}} to use mathematical functions not implemented by the package.
}


\value{
SpatRaster or SpatExtent
}


\examples{
r1 <- rast(ncols=10, nrows=10)
v <- runif(ncell(r1))
v[10:20] <- NA
values(r1) <- v
r2 <- rast(r1)
values(r2) <- 1:ncell(r2) / ncell(r2)
r3 <- r1 + r2
r2 <- r1 / 10
r3 <- r1 * (r2 - 1 / r2)

b <- c(r1, r2, r3)
b2 <- b * 10

### SpatExtent methods
x <- ext(0.1, 2.2, 0, 3)
y <- ext(-2, 1, -2,2)
# union
x + y
# intersection
x * y

e <- x %% 2
e
e * 2
e / 2
e + 1
e - 1
}

\keyword{methods}
\keyword{math}
\keyword{spatial}

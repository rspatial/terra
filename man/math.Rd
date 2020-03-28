\name{math}

\docType{methods}

\alias{Arith-methods}
\alias{Compare-methods}
\alias{Math-methods}
\alias{Math2-methods}

\alias{Arith,SpatRaster,SpatRaster-method}
\alias{Arith,SpatRaster,missing-method}
\alias{Arith,numeric,SpatRaster-method}
\alias{Arith,SpatRaster,numeric-method}
\alias{Arith,missing,SpatRaster-method}
\alias{Compare,SpatRaster,SpatRaster-method}
\alias{Compare,numeric,SpatRaster-method}
\alias{Compare,SpatRaster,numeric-method}
\alias{Compare,SpatRaster,math-method}
\alias{Math,SpatRaster-method}
\alias{Math2,SpatRaster-method}
\alias{Math,SpatExtent-method}
\alias{Math2,SpatExtent-method}
\alias{Compare,SpatExtent,SpatExtent-method}

\alias{is.na}
\alias{is.na,SpatRaster-method}

\title{Mathematical operations with SpatRaster objects}

\description{

Standard arithmetic operators for computations with SpatRaster objects and numeric values. The following operators and methods are available for SpatRaster:

\bold{Arith}: \code{ +, -, *, /, ^, \%\%, \%/\% }

\bold{Compare}: \code{==, !=, >, <,  <=, >=, is.na}

\bold{\code{\link[base]{Summary}}}: \code{"max", "min", "range", "prod", "sum", "any", "all"}
 
\bold{Math}: \code{"abs", "sign", "sqrt", "ceiling", "floor", "trunc",  "cummax", "cummin", "cumprod", "cumsum", "log", "log10", "log2", "log1p", "acos", "acosh", "asin", "asinh", "atan", "atanh", "exp", "expm1", "cos", "cosh", "sin", "sinh", "tan", "tanh"}
 
\bold{Math2}: \code{"round", "signif"}

If multiple SpatRaster objects are used, these must have the same extent and resolution. Operators are applied on a cell by cell basis.

For \bold{SpatExtent} the following methods have been implemented: \code{"round", "floor", "ceil", "=="}

}


\usage{
\S4method{is.na}{SpatRaster}(x)
}


\section{Note}{
You can use the somewhat more flexible \code{\link{app}} method instead of the Math-methods.
}


\value{
SpatRaster
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
s <- sqrt(b2)
round(s, 1)

x <- is.na(s)
}

\keyword{methods}
\keyword{math}
\keyword{spatial}

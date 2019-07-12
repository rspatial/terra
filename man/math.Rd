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

\title{Mathematical operations with SpatRaster objects}

\description{
Standard arithmetic operators for computations with SpatRaster objects and numeric values. The following operators and methods are available:

\bold{Arith}: \code{ +, -, *, /, ^, \%\%, \%/\% }

\bold{Compare}: \code{==, !=, >, <,  <=, >=}

\bold{\code{\link{Summary}}}: \code{"max", "min", "range", "prod", "sum", "any", "all"}
 
\bold{Math}: \code{"abs", "sign", "sqrt", "ceiling", "floor", "trunc",  "cummax", "cummin", "cumprod", "cumsum", "log", "log10", "log2", "log1p", "acos", "acosh", "asin", "asinh", "atan", "atanh", "exp", "expm1", "cos", "cosh", "sin", "sinh", "tan", "tanh"}
 
\bold{Math2}: \code{"round", "signif"}
 
If multiple SpatRaster objects are used, these must have the same extent and resolution. Operators are applied on a cell by cell basis.
}


\section{Note}{
You can use the somewhat more flexible \code{\link{reduce}} method instead of the Math-methods.
}


\value{
SpatRaster
}

\examples{
r1 <- rast(ncols=10, nrows=10)
values(r1) <- runif(ncell(r1))
r2 <- rast(r1)
values(r2) <- 1:ncell(r2) / ncell(r2)
r3 <- r1 + r2
r2 <- r1 / 10
r3 <- r1 * (r2 - 1 / r2)


b <- c(r1, r2, r3)
b2 <- b * 10
s <- sqrt(b2)
round(s, 1)
}

\keyword{methods}
\keyword{math}
\keyword{spatial}

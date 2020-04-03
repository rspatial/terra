\name{math}

\docType{methods}

\alias{Arith-methods}
\alias{Compare-methods}
\alias{Logic-methods}
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
\alias{Logic,SpatRaster,SpatRaster-method}
\alias{Logic,SpatRaster,numeric-method}
\alias{Logic,numeric,SpatRaster-method}
\alias{Math,SpatRaster-method}
\alias{Math2,SpatRaster-method}
\alias{Math,SpatExtent-method}
\alias{Math2,SpatExtent-method}
\alias{Compare,SpatExtent,SpatExtent-method}

\alias{is.na,SpatRaster-method}
\alias{is.nan,SpatRaster-method}
\alias{!,SpatRaster-method}
\alias{is.finite,SpatRaster-method}
\alias{is.infinite,SpatRaster-method}


\title{Arithmetic, logical and general mathematical methods}

\description{

Standard operators and mathematical methods for computations with SpatRaster objects. Computations are local (applied on a cell by cell basis). If multiple SpatRaster objects are used, these must have the same extent and resolution. These have been implemented:

\bold{Arith}: \code{ +, -, *, /, ^, \%\%, \%/\% }

\bold{Compare}: \code{ ==, !=, >, <,  <=, >=, is.na, is.nan, is.finite, is.infinite}

The terra package does not distunguish between \code{NA} (not available) and \code{NaN} (not a number) and that both are represented as \code{NaN}.  

\bold{Logical}: \code{!, &, |}

\bold{Summary}: \code{"max", "min", "range", "prod", "sum", "any", "all"}
 
\bold{Math}: \code{"abs", "sign", "sqrt", "ceiling", "floor", "trunc",  "cummax", "cummin", "cumprod", "cumsum", "log", "log10", "log2", "log1p", "acos", "acosh", "asin", "asinh", "atan", "atanh", "exp", "expm1", "cos", "cosh", "sin", "sinh", "tan", "tanh", "round", "signif"}


For \bold{SpatExtent} the following methods have been implemented: \code{"round", "floor", "ceil", "=="}
}


\usage{
\S4method{is.na}{SpatRaster}(x)
}

\arguments{
  \item{x}{SpatRaster}
}

\section{seealso}{
\code{\link{ifel}} to convieniently combine operations and \code{\link{app}} to use mathematical functions not implemented by the package.
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

max(s)
max(s, na.rm=TRUE)

x <- is.na(s)

y <- which.max(s)
}

\keyword{methods}
\keyword{math}
\keyword{spatial}

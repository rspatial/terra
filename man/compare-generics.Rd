\name{Compare-methods}

\docType{methods}

\alias{Compare-methods}
\alias{Logic-methods}
\alias{logic}
\alias{compare}

\alias{Compare,SpatRaster,SpatRaster-method}
\alias{Compare,numeric,SpatRaster-method}
\alias{Compare,SpatRaster,numeric-method}
\alias{Compare,SpatRaster,character-method}
\alias{Compare,SpatRaster,matrix-method}
\alias{Compare,matrix,SpatRaster-method}


\alias{Logic,SpatRaster,SpatRaster-method}
\alias{Logic,SpatRaster,numeric-method}
\alias{Logic,numeric,SpatRaster-method}
\alias{Logic,SpatRaster,logical-method}
\alias{Logic,logical,SpatRaster-method}

\alias{Compare,SpatExtent,SpatExtent-method}

\alias{is.na,SpatRaster-method}
\alias{is.nan,SpatRaster-method}
\alias{!,SpatRaster-method}
\alias{is.finite,SpatRaster-method}
\alias{is.infinite,SpatRaster-method}

\alias{logic,SpatRaster-method}
\alias{compare,SpatRaster-method}

\title{Compare and logical methods}

\description{

Standard comparison and logical operators for computations with SpatRasters. Computations are local (applied on a cell by cell basis). If multiple SpatRaster objects are used, these must have the same geometry (extent and resolution). These operators have been implemented:

\bold{Logical}: \code{!, &, |, isTRUE, isFALSE}

\bold{Compare}: \code{ ==, !=, >, <,  <=, >=, is.na, is.nan, is.finite, is.infinite}

See \code{\link{not.na}} for the inverse of \code{is.na}, and \code{\link{noNA}} to detect cells with missing value across layers.

The \code{compare} and \code{logic} methods implement these operators in a method that can return \code{NA} istead of \code{FALSE} and allows for setting an output filename.

The terra package does not distinguish between \code{NA} (not available) and \code{NaN} (not a number). In most cases this state is represented by \code{NaN}.  

If you use a SpatRaster with a vector of multiple numbers, each element in the vector is considered a layer (with a constant value). If you use a SpatRaster with a matrix, the number of columns of the matrix must match the number of layers of the SpatRaster. The rows are used to match the cells. That is, if there are two rows, these match cells 1 and 2, and they are recycled to 3 and 4, etc. 

The following method has been implemented for \bold{(SpatExtent, SpatExtent)}: \code{==}
}


\usage{
\S4method{compare}{SpatRaster}(x, y, oper, falseNA=FALSE, filename="", overwrite=FALSE, ...)

\S4method{logic}{SpatRaster}(x, oper, falseNA=FALSE, filename="", overwrite=FALSE, ...)
}

\arguments{
\item{x}{SpatRaster}
\item{y}{SpatRaster or numeric}
\item{oper}{character. Operator name. For \code{compare} this can be one of \code{"==", "!=", ">", "<", ">=", "<="} and for \code{logic} it can be one of \code{"!", "is.na", "allNA", "noNA", "is.infinite", "is.finite", "iSTRUE", "isFALSE"}}
\item{falseNA}{logical. Should the result be \code{TRUE, NA} instead of \code{TRUE, FALSE}?}
\item{filename}{character. Output filename}
\item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
\item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}


\seealso{
\code{\link{all.equal}}, \link{Arith-methods}. See \code{\link{ifel}} to conveniently combine operations and \code{\link{Math-methods}} or \code{\link{app}} to apply any R function to a SpatRaster. 
}


\value{
SpatRaster or SpatExtent
}


\examples{
r1 <- rast(ncols=10, nrows=10)
values(r1) <- runif(ncell(r1))
r1[10:20] <- NA
r2 <- rast(r1)
values(r2) <- 1:ncell(r2) / ncell(r2)

x <- is.na(r1)
!x
r1 == r2 
compare(r1, r2, "==")
compare(r1, r2, "==", TRUE)
}

\keyword{methods}
\keyword{math}
\keyword{spatial}

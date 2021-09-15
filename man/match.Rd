\name{match}

\docType{methods}

\alias{match}
\alias{match,SpatRaster-method}
\alias{\%in\%}
\alias{\%in\%,SpatRaster-method}

\title{Value matching for SpatRasters}

\description{
\code{match} returns a SpatRaster with the position of the matched values. The cell values are the index of the table argument. 

\code{\%in\%} returns a 0/1 (FALSE/TRUE) SpatRaster indicating if the cells values were matched or not.

}

\usage{
match(x, table, nomatch = NA_integer_, incomparables = NULL)

x \%in\% table
}

\arguments{
 \item{x}{SpatRaster}
 \item{table}{vector of the values to be matched against}
 \item{nomatch}{the value to be returned in the case when no match is found. Note that it is coerced to integer}
 \item{incomparables}{a vector of values that cannot be matched. Any value in x matching a value in this vector is assigned the nomatch value. For historical reasons, FALSE is equivalent to NULL}
}

\value{
SpatRaster
}

\seealso{ \code{\link{app}, \link[base]{match}} }




\examples{
r <- rast(nrows=10, ncols=10)
values(r) <- 1:100
m <- match(r, c(5:10, 50:55))
n <- r \%in\% c(5:10, 50:55)
}

\keyword{spatial}
\keyword{methods}

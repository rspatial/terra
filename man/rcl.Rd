\name{rcl}

\alias{rcl}
\alias{rcl,SpatRaster-method}


\title{Combine row, column, and layer numbers}

\description{
Get a matrix with the combination of row, column, and layer numbers
}

\usage{
\S4method{rcl}{SpatRaster}(x, row=NULL, col=NULL, lyr=NULL)
}

\arguments{
  \item{x}{SpatRaster}
  \item{row}{positive integer that are row number(s), a list thereof, or NULL for all rows}
  \item{col}{as above for columns}
  \item{lyr}{as above for layers}
}

\details{
If a list is used for at least one of \code{row}, \code{col} or \code{lyr}, these are evaluated in parallel. That is combinations are made for each list element, not across list elements. If, in this case another argument is not a list it has to have either length 1 (used for all cases) or have the same length as the (longest) list, in which case the value is coerced into a list with \code{as.list}

If multiple arguments are a list but they have different lengths, theyare recycled to the longest list. 
}

\value{
matrix
}


\seealso{
\code{\link{rowColCombine}}, \code{\link{cellFromRowCol}}
}



\examples{
x <- rast(ncol=5, nrow=5, nlyr=2)
values(x) <- 1:size(x)

rcl(x, 1, 2:3, 1:2)

i <- rcl(x, 1, list(1:2, 3:4), 1:2)
i

# get the values for these cells
x[i]
}

\keyword{spatial}

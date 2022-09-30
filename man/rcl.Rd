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
  \item{row}{positive integer that are row number(s) or NULL for all rows}
  \item{col}{as above for columns}
  \item{col}{as above for layers}
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

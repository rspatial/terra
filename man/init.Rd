\name{initialize}

\alias{init}
\alias{init,SpatRaster-method}

\title{Initialize a SpatRaster with values}

\description{
Create a SpatRaster with values reflecting a cell property: "x", "y", "col", "row", "cell" or "chess". Alternatively, a function can be used. In that case, cell values are initialized without reference to pre-existing values. E.g., initialize with a random number (\code{fun=\link{runif}}). While there are more direct ways of achieving this for small objects (see examples) for which a vector with all values can be created in memory, the \code{init} function will also work for SpatRasters with many cells.
}

\usage{
\S4method{init}{SpatRaster}(x, fun, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{fun}{function to be applied. This must be a either single number, multiple numbers, a function, or one of a set of known character values. A function must take the number of cells as a single argument to return a vector of values with a length equal to the number of cells, such as \code{fun=runif}. Allowed character values are "x", "y", "row", "col", "cell", and "chess" to get the x or y coordinate, row, col or cell number or a chessboard pattern (alternating 0 and 1 values)}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}


\value{
SpatRaster
}

\examples{
r <- rast(ncols=10, nrows=5, xmin=0, xmax=10, ymin=0, ymax=5)
x <- init(r, fun="cell")
y <- init(r, fun=runif)

# initialize with a single value 
z <- init(r, fun=8)
}

\keyword{spatial}


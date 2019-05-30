\name{initialize}

\alias{init}

\title{Initialize a SpatRaster with values}

\description{
Create a SpatRaster with values reflecting a cell property: 'x', 'y', 'col', 'row', 'cell' or 'chess'. Alternatively, a function can be used. In that case, cell values are initialized without reference to pre-existing values. E.g., initialize with a random number (\code{fun=\link{runif}}). While there are more direct ways of achieving this for small objects (see examples) for which a vector with all values can be created in memory, the \code{init} function will also work for Raster* objects with many cells.
}

\usage{
\S4method{init}{Raster}(x, fun, filename="", ...)
}

\arguments{
  \item{x}{Raster* object}
  \item{fun}{function to be applied. This must be a function that can take the number of cells as a single argument to return a vector of values with a length equal to the number of cells, such as \code{fun=runif}. You can also supply one of the following character values: 'x', 'y', 'row', 'col', or 'cell' to get the x or coordinate, row, col or cell number; you can also use 'chess', to get a chessboard pattern}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{Additional arguments. None implemented}  
}


\value{
SpatRaster
}

\examples{
r <- rast(ncol=10, nrow=5, xmn=0, xmx=10, ymn=0, ymx=5)
x <- init(r, fun="cell")
y <- init(r, fun=runif)

# initialize with a single value 
z <- init(r, fun=function(n){rep(3, n)})

}

\keyword{spatial}


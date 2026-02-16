\name{chunk}

\alias{chunk}
\alias{chunk,SpatRaster-method}

\title{Make a SpatRaster method memory-safe}

\description{
This method allows for running a function that takes a SpatRaster as first argument in chunks (tiles). This can be useful if the functions is not memory-safe, typically because it reads all the raster cell values into memory. 

This method is not designed to be especially efficient, and there might be more efficient ways to accomplish what the goal of the function that is not memory-safe.

Also, some functions must have access to all cells at once to be valid. In those cases, chunk would return incorrect results.
}

\usage{
\S4method{chunk}{SpatRaster}(x, fun, ..., n=NULL, buffer=0, filename="", wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{fun}{function that takes a SpatRaster as first argument}
  \item{...}{additional arguments for fun}
  \item{n}{NULL or positive integer to specifying the number of rows and columns for each chunk (or 2 numbers for a different number of rows and columns, as in \code{\link{getTileExtents}}}
  \item{buffer}{integer. The number of additional rows and columns added to each tile. Can be a single number, or two numbers to specify a separate number of rows and columns. This allows for creating overlapping tiles that can be used for computing spatial context dependent values with e.g. \code{\link{focal}}. The expansion is only inside \code{x}, no rows or columns outside of \code{x} are added}
  \item{filename}{character. Output filename}
  \item{wopt}{list with additional arguments for writing files as in \code{\link{writeRaster}}}
}


\value{SpatRaster }


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))
f <- function(x, a = 0) {
	print("chunk")
	sum(x) + a
}

x <- chunk(s, f, a=100)

}

\keyword{spatial}


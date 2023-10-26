\name{sieve}

\alias{sieve}
\alias{sieve,SpatRaster-method}

\title{Sieve filter}

\description{
Apply a sieve filter. That is, remove "noise", by changing small clumps of cells with a value that is different from the surrounding cells, to the value of the largest neighboring clump.

Note that the input raster data is rounded to integers.
}

\usage{
\S4method{sieve}{SpatRaster}(x, threshold, directions=8, filename="", ...) 
}

\arguments{
  \item{x}{SpatRaster, single layer with elevation values. Values should have the same unit as the map units}
  
  \item{threshold}{positive integer. Only clumps smaller than this threshold will be removed}
  
  \item{directions}{numeric to indicate which cells are connected. Either \code{4} to only consider the horizontal and vertical neighbors ("rook"), or \code{8} to consider the vertical, horizontal and diagonal neighbors}

  \item{filename}{character. Output filename}
  \item{...}{Options for writing files as in \code{\link{writeRaster}}}
}


\seealso{\code{\link{focal}}}


\examples{
r <- rast(nrows=18, ncols=18, xmin=0, vals=0, crs="local")
r[2, 5] <- 1
r[5:8, 2:3] <- 2
r[7:12, 10:15] <- 3
r[15:16, 15:18] <- 4
freq(r, bylayer=FALSE)

x <- sieve(r, 8)
y <- sieve(r, 9)
}

\keyword{spatial}

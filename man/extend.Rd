\name{extend}

\alias{extend}
\alias{extend,SpatRaster-method}
\alias{extend,SpatExtent-method}


\title{Extend}

\description{
Enlarge the spatial extent of a SpatRaster. See \code{\link{crop}} if you (also) want to remove rows or columns. 

Note that you can only enlarge SpatRasters with entire rows and columns. Therefore, the extent of the output SpatRaster may not be exactly the same as the requested. Depending on argument \code{snap} it may be a bit smaller or larger. 

You can also enlarge a SpatExtent with this method, or with an algebraic notation (see examples)
}

\usage{
\S4method{extend}{SpatRaster}(x, y, snap="near", fill=NA, filename="", overwrite=FALSE, ...) 

\S4method{extend}{SpatExtent}(x, y) 
}

\arguments{
  \item{x}{SpatRaster or SpatExtent}
  \item{y}{If \code{x} is a SpatRaster, \code{y} should be a SpatExtent, or an object from which it can be extracted (such as SpatRaster and SpatVector objects). Alternatively, you can provide one, two or four non-negative integers indicating the number of rows and columns that need to be added at each side (a single positive integer when the number of rows and columns to be added is equal; or 2 number (columns, rows), or four (left column, right column, bottom row, top row).
  If \code{x} is a SpatExtent, \code{y} should likewise be a numeric vector of 1, 2, or 4 elements}
  
  \item{snap}{character. One of "near", "in", or "out". Used to align \code{y} to the geometry of \code{x}} 
  \item{fill}{numeric. The value used to for the new raster cells}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster or SpatExtent
}

\seealso{\code{\link{crop}}, \code{\link{merge}}, \code{\link{ext}}, \code{\link{resample}}}

\examples{
r <- rast(xmin=-150, xmax=-120, ymin=30, ymax=60, ncols=36, nrows=18)
values(r) <- 1:ncell(r)
e <- ext(-180, -100, 40, 70)
re <- extend(r, e)

# extend with a number of rows and columns (at each side)
re2 <- extend(r, c(2,10))

# SpatExtent
e <- ext(r)
e
extend(e, 10)
extend(e, c(10, -10, 0, 20))


# add 10 columns / rows on all sides
e + 10
# double extent
e * 2
# increase extent by 25%
e * 1.25
}

\keyword{spatial}

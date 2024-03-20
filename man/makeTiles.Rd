\name{makeTiles}

\docType{methods}


\alias{getTileExtents}
\alias{getTileExtents,SpatRaster-method}

\alias{makeTiles}
\alias{makeTiles,SpatRaster-method}


\title{Make tiles or get their extents}

\description{ 
Divide a SpatRaster into "tiles". The cells of another SpatRaster (normally with a much lower resolution) or a SpatVector with polygon geometry can be used to define the tiles. You can also provide one or two numbers to indicate the number of rows and columns per tile. 

\code{getTileExtents} returns the extents of the (virtual) tiles, while \code{makeTiles} creates files for the tiles and returns their filenames. 
}

\usage{
\S4method{makeTiles}{SpatRaster}(x, y, filename="tile_.tif", extend=FALSE,
		na.rm=FALSE, buffer=0, overwrite=FALSE, ...)

\S4method{getTileExtents}{SpatRaster}(x, y, extend=FALSE, buffer=0)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatRaster or SpatVector defining the zones; or numeric specifying the number of rows and columns for each zone (1 or 2 numbers if the number of rows and columns is not the same)}
  \item{filename}{character. Output filename template. Filenames will be altered by adding the tile number for each tile}
  \item{extend}{logical. If \code{TRUE}, the extent of \code{y} is expanded to assure that it covers all of \code{x}}
  \item{na.rm}{logical. If \code{TRUE}, tiles with only missing values are ignored}
  \item{buffer}{integer. The number of additional rows and columns added to each tile. Can be a single number, or two numbers to specify a separate number of rows and columns. This allows for creating overlapping tiles that can be used for computing spatial context dependent values with e.g. \code{\link{focal}}. The expansion is only inside \code{x}, no rows or columns outside of \code{x} are added}
  \item{overwrite}{logical. If \code{TRUE}, existing tiles are overwritten; otherwise they are skipped (without error or warning)}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
character (filenames) or matrix (extents)
}


\seealso{
\code{\link{vrt}} to create a virtual raster from tiles and \code{\link{crop}} for sub-setting arbitrary parts of a SpatRaster.
}



\examples{
r <- rast(ncols=100, nrows=100)
values(r) <- 1:ncell(r)
x <- rast(ncols=2, nrows=2)

getTileExtents(r, x)
getTileExtents(r, x, buffer=3)


filename <- paste0(tempfile(), "_.tif")
ff <- makeTiles(r, x, filename)
ff

vrt(ff)
}


\keyword{methods}
\keyword{spatial}

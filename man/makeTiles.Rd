\name{makeTiles}

\docType{methods}

\alias{tiles}
\alias{tiles,SpatRaster-method}
\alias{makeTiles}
\alias{makeTiles,SpatRaster-method}


\title{Make tiles}

\description{ 
Divide a SpatRaster into "tiles". The cells of another SpatRaster (normally with a much lower resolution) or a SpatVector with polygon geometry can be used to define the tiles. You can also provide one or two numbers to indicate the number of rows and columns per tile. 
}

\usage{
\S4method{makeTiles}{SpatRaster}(x, y, filename="tile_.tif", extend=FALSE,
		na.rm=FALSE, buffer=0, overwrite=FALSE, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatRaster, SpatVector, or numeric}
  \item{filename}{character. Output filename template. Filenames will be altered by adding the tile number for each tile}
  \item{extend}{logical. If \code{TRUE}, the extent of \code{y} is expanded to assure that it covers all of \code{x}}
  \item{na.rm}{logical. If \code{TRUE}, tiles with only missing values are ignored}
  \item{buffer}{integer. The number of cells to add around each tile. Ignored if \code{y} is a SpatVector. Can be a single number, or two numbers to specify a separate buffer for the rows and for the columns. This allows for creating overlapping tiles that can be used for computing spatial context dependent values with e.g. \code{\link{focal}}. The expansion is only inside \code{x}, no rows or columns outside of \code{x} are added}
  \item{overwrite}{logical. If \code{TRUE}, existing tiles are overwritten; otherwise they are skipped (without error or warning)}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
character (filenames)
}


\seealso{
\code{\link{vrt}} to create a virtual raster from tiles
}



\examples{
r <- rast(ncols=100, nrows=100)
values(r) <- 1:ncell(r)
x <- rast(ncols=2, nrows=2)
filename <- paste0(tempfile(), "_.tif")
ff <- makeTiles(r, x, filename)
ff

vrt(ff)
}


\keyword{methods}
\keyword{spatial}

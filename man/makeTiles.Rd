\name{makeTiles}

\docType{methods}

\alias{tiles}
\alias{tiles,SpatRaster-method}
\alias{makeTiles}
\alias{makeTiles,SpatRaster-method}


\title{Make tiles}

\description{ 
Divide a SpatRaster into "tiles". The cell of another SpatRaster (normally with a much lower resolution) are used to define the tiles.
}

\usage{
\S4method{makeTiles}{SpatRaster}(x, y, expand=FALSE, filename="tile_.tif", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatRaster or SpatVector}
  \item{expand}{logical. If \code{TRUE}, the extent of \code{y} is expanded to assure that it covers all of \code{x}}
  \item{filename}{character. Output filename template. Filenames will be altered by adding the tile number for each tile}
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

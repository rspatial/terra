\name{vrt}

\docType{methods}

\alias{vrt}
\alias{vrt,character-method}


\title{Virtual Raster Tiles}

\description{ 
Create a Virtual Raster Tiles (VRT) dataset from a collection of file-based raster datasets. 
}

\usage{
\S4method{vrt}{character}(x, filename="", overwrite=FALSE)
}

\arguments{
  \item{x}{character. Filenames of raster "tiles". See \code{\link{tiles}}}
  \item{filename}{character. Output VRT filename}
  \item{overwrite}{logical. Should \code{filename} be overwritten if it exists?}
}

\value{
SpatRaster
}

\seealso{\code{\link{makeVRT}} to create a .vrt file for a file without a header}

\examples{
r <- rast(ncols=100, nrows=100)
values(r) <- 1:ncell(r)
x <- rast(ncols=2, nrows=2)
filename <- paste0(tempfile(), "_.tif")
ff <- makeTiles(r, x, filename)
ff

vrtfile <- paste0(tempfile(), ".vrt")
v <- vrt(ff, vrtfile)
head(readLines(vrtfile))
v
}


\keyword{methods}
\keyword{spatial}

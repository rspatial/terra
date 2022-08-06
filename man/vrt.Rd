\name{vrt}

\docType{methods}

\alias{vrt}
\alias{vrt,character-method}


\title{Virtual Raster Tiles}

\description{ 
Create a Virtual Raster Tiles (VRT) dataset from a collection of file-based raster datasets. 
}

\usage{
\S4method{vrt}{character}(x, filename="", options=NULL, overwrite=FALSE)
}

\arguments{
  \item{x}{character. Filenames of raster "tiles". See \code{\link{tiles}}}
  \item{filename}{character. Output VRT filename}  
  \item{options}{character. All arguments as separate vector elements. Options as for \href{https://gdal.org/programs/gdalbuildvrt.html}{gdalbuildvrt}}
  \item{overwrite}{logical. Should \code{filename} be overwritten if it exists?}
}

\value{
SpatRaster
}

\seealso{
\code{\link{makeTiles}} to create tiles; \code{\link{makeVRT}} to create a .vrt file for a file without a header
}

\examples{
r <- rast(ncols=100, nrows=100)
values(r) <- 1:ncell(r)
x <- rast(ncols=2, nrows=2)
filename <- paste0(tempfile(), "_.tif")
ff <- makeTiles(r, x, filename)
ff

vrtfile <- paste0(tempfile(), ".vrt")
# output in lower resolution
v <- vrt(ff, vrtfile, options = c("-tr", 5, 5))
head(readLines(vrtfile))
v
}


\keyword{methods}
\keyword{spatial}

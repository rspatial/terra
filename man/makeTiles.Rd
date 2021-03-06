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
\S4method{makeTiles}{SpatRaster}(x, y, filename="tile_.tif", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatRaster or SpatVector}
  \item{filename}{character. Output filename template. Filenames will be altered by adding the tilenumber for each tile}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
character (filenames)
}


\examples{
r <- rast(ncols=100, nrows=100)
values(r) <- 1:ncell(r)
x <- rast(ncols=2, nrows=2)
filename <- paste0(tempfile(), "_.tif")
ff <- makeTiles(r, x, filename)
}


\keyword{methods}
\keyword{spatial}

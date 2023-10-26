\name{wrapCache}

\alias{wrapCache}

\alias{wrapCache,SpatRaster-method}


\title{SpatRaster wrap with caching options}

\description{
Use \code{wrap} to pack a SpatRaster with caching options. See \code{\link{wrap}} for the general approach that is easier and better to use in most cases. 

This method allows for specifying a folder, or filenames, to cache all sources of a SpatRaster in a specific location (on disk).
}

\usage{
\S4method{wrapCache}{SpatRaster}(x, filename=NULL, path=NULL, overwrite=FALSE, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{filename}{character. A single filename, or one filename per SpatRaster data source. If not \code{NULL}, the raster sources are saved in these files}
  \item{path}{character. If not \code{NULL}, the path where raster sources will be saved. Ignored if filenames is not \code{NULL}}
  \item{overwrite}{Should existing files be overwritten when \code{files} or \code{path} is not \code{NULL}? If this value is not \code{TRUE} or \code{FALSE}, only files that do not exist are created}
  \item{...}{Additional arguments for \code{writeRaster}. Only used for raster sources that are in memory, as other sources are cached by copying the files}  
}

\value{
PackedSpatRaster
}

\seealso{\code{\link{wrap}}, \code{\link{unwrap}}}

\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)

x <- wrapCache(r, path=tempdir())
x

}


\keyword{ spatial }
\keyword{ methods }


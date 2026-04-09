\name{update}

\alias{update}
\alias{update,SpatRaster-method}

  
\title{Change values in a file}

\description{
Change the contents of a file that is the data source of a SpatRaster. BE CAREFUL as you are overwriting values in an existing file.
}

\usage{
\S4method{update}{SpatRaster}(object, crs=FALSE, extent=FALSE, names=FALSE)
}

\arguments{
  \item{object}{SpatRaster}
  \item{crs}{logical. Should the coordinate reference system be updated?}
  \item{extent}{logical. Should the extent be updated?}
  \item{names}{logical. Should the names be updated?}
}

\value{SpatRaster (invisibly)}


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
fname <- paste0(tempfile(), ".tif")
x <- writeRaster(s*1, fname)

ext(x) <- ext(x) + 1
crs(x) <- "+proj=utm +zone=1"
names(x) <- LETTERS[1:3]"

update(x, crs=TRUE, extent=TRUE)

rast(fname)
}

\keyword{spatial}


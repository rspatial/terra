\name{update}

\alias{update}
\alias{update,SpatRaster-method}

  
\title{Change values in a file}

\description{
Change the contents of a file that is the data source of a SpatRaster. BE CAREFUL as you are overwriting values in an existing file.
}

\usage{
\S4method{update}{SpatRaster}(object, names=FALSE, crs=FALSE, extent=FALSE)
}

\arguments{
  \item{object}{SpatRaster}
  \item{names}{logical. Should the names be updated?}
  \item{crs}{logical. Should the coordinate reference system be updated?}
  \item{extent}{logical. Should the extent be updated?}
}

\value{SpatRaster (invisibly)}


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
fname <- paste0(tempfile(), ".tif")
x <- writeRaster(s, fname)
names(x) <- c("A", "B", "C")
ext(x) <- ext(x) + 1
crs(x) <- "+proj=utm +zone=1"

update(x, names=TRUE, crs=TRUE, extent=TRUE)

rast(fname)
}

\keyword{spatial}


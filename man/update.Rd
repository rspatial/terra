\name{update}

\alias{update}
\alias{update,SpatRaster-method}

  
\title{Change values in a file}

\description{
Change the contents of a file that is the data source of a SpatRaster. BE CAREFUL as you are overwriting values in an existing file.
}

\usage{
\S4method{update}{SpatRaster}(object, names=FALSE)
}

\arguments{
  \item{object}{SpatRaster}
  \item{names}{logical. Should the names be updated?}
}

\value{SpatRaster (invisibly)}


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
fname <- paste0(tempfile(), ".tif")
x <- writeRaster(s, fname)
names(x) <- c("A", "B", "C")
update(x, names=TRUE)

r <- rast(fname)
names(r)
}

\keyword{spatial}


\name{update}

\alias{update.names}
\alias{update.names,SpatRaster-method}

  
\title{Change values in a file}

\description{
The "update." replacement methods change the contents of an existing file. BE CAREFUL. 

\code{update.names} writes the (changed) names of SpatRaster \code(x) to the file of the data source.
}

\usage{
\S4method{update.names}{SpatRaster}()
}

\arguments{
  \item{x}{SpatRaster}
}

\value{SpatRaster (invisibly)}


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
fname <- paste0(tempfile(), ".tif")
x <- writeRaster(s, fname)
names(x) <- c("A", "B", "C")
update.names(x)

r <- rast(fname)
names(r)
}

\keyword{spatial}


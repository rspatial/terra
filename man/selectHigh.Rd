\name{selectHighest}

\docType{methods}

\alias{selectHighest}
\alias{selectHighest,SpatRaster-method}

\title{select cells with high or low values}

\description{
Identify n cells that have the highest of lowest values in the first layer of a SpatRaster. 
}

\usage{
\S4method{selectHighest}{SpatRaster}(x, n, low=FALSE)
}

\arguments{
  \item{x}{SpatRaster. Only the first layer is processed}
  \item{n}{The number of cells to select}
  \item{low}{logical. If \code{TRUE}, the lowest values are selected instead of the highest values}
}


\value{
SpatRaster
}

\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
x <- selectHighest(r, 1000)
y <- selectHighest(r, 1000, TRUE)

m <- merge(y-1, x)
levels(m) <- data.frame(id=0:1, elevation=c("low", "high"))
plot(m)
}

\keyword{spatial}

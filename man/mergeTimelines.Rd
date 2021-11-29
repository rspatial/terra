\name{mergeTimelines}

\docType{methods}

\alias{mergeTimelines}
\alias{mergeTimelines,SpatRasterCollection-method}

\title{ 
merge SpatRasters by timelines
}

\description{
merge SpatRasters by timelines
}

\usage{
\S4method{mergeTimelines}{SpatRasterCollection}(x, fun=mean, filename="", ...)
}

\arguments{
  \item{x}{SpatRasterDataset}
  \item{fun}{reducing function}
  \item{filename}{character. Output filename}
  \item{...}{list with named options for writing files as in \code{\link{writeRaster}}}
} 

\value{
SpatRaster
}


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
s <- c(s, s)
time(s) <- as.Date("2001-01-01") + 0:5
s1 <- s/10
time(s1) <- as.Date("2001-01-07") + 0:5
s2 <- s*10
time(s2) <- as.Date("2001-01-05") + 0:5
x <- sds(s, s1, s2)

m <- mergeTimelines(x, mean)
}

\keyword{methods}
\keyword{spatial}


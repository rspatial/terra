\name{mergeTimelines}

\docType{methods}

\alias{mergeTimelines}
\alias{mergeTimelines,SpatRasterDataset-method}

\title{ 
merge SpatRasters by timelines to create a single timeseries
}

\description{
merge SpatRasters with partly overlapping timestamps to create a single time series. If there is no overlap there is no point in using this functoin (use \code{\link{c}} instead).

Time gaps are not filled 
}

\usage{
\S4method{mergeTimelines}{SpatRasterDataset}(x, fun=mean, filename="", ...)
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
r <- rast(system.file("ex/logo.tif", package="terra"))   
s1 <- c(r, r)
time(s1) <- as.Date("2001-01-01") + 0:5
s1 <- s1/10
time(s1) <- as.Date("2001-01-07") + 0:5
s2 <- s1*10
time(s2) <- as.Date("2001-01-05") + 0:5
x <- sds(s1, s1, s2)

m <- mergeTimelines(x, mean)
}

\keyword{methods}
\keyword{spatial}


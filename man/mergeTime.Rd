\name{mergeTime}

\docType{methods}

\alias{mergeTime}
\alias{mergeTime,SpatRasterDataset-method}

\title{ 
merge SpatRasters by timelines to create a single timeseries
}

\description{
Combine SpatRasters with partly overlapping time-stamps to create a single time series. If there is no overlap between the SpatRasters there is no point in using this function (use \code{\link{c}} instead). 

Also note that time gaps are not filled. You can use \code{\link{fillTime}} to do that.
}

\usage{
\S4method{mergeTime}{SpatRasterDataset}(x, fun=mean, filename="", ...)
}

\arguments{
  \item{x}{SpatRasterDataset}
  \item{fun}{A function that reduces a vector to a single number, such as \code{mean} or \code{min}}
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

m <- mergeTime(x, mean)
}

\keyword{methods}
\keyword{spatial}


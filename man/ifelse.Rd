\name{ifel}

\alias{ifel}
\alias{ifel,SpatRaster-method}

\title{ifelse for SpatRasters}

\description{
Implementation of \code{\link[base]{ifelse}} for SpatRasters. This method allows for a concise expression of what can otherwise be achieved with a combination of \code{\link{classify}}, \code{\link{mask}}, and \code{\link{cover}}. 

\code{ifel} is an \code{R} equivalent to the \code{Con} method in ArcGIS (arcpy). 
}

\usage{
\S4method{ifel}{SpatRaster}(test, yes, no, filename="", ...)
}

\arguments{
  \item{test}{SpatRaster with logical (TRUE/FALSE) values}
  \item{yes}{SpatRaster or numeric}
  \item{no}{SpatRaster or numeric}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\examples{
r <- rast(nrows=5, ncols=5, xmin=0, xmax=1, ymin=0, ymax=1)
values(r) <- c(-10:0, NA, NA, NA, 0:10)

x <- ifel(r > 1, 1, r)
# same as 
a <- classify(r, cbind(1, Inf, 1))
# or
b <- app(r, fun=function(i) {i[i > 1] <- 1; i})
# or 
d <- clamp(r, -Inf, 1)
# or (not recommended for large datasets)
e <- r
e[e>1] <- 1

## other examples
f <- ifel(is.na(r), 100, r)

z <- ifel(r > -2 & r < 2, 100, 0)

# nested expressions
y <- ifel(r > 1, 1, ifel(r < -1, -1, r))

k <- ifel(r > 0, r+10, ifel(r < 0, r-10, 3))
}

\keyword{spatial}

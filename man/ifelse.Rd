\name{ifel}

\alias{ifel}
\alias{ifel,SpatRaster-method}

\title{ifelse for SpatRasters}

\description{
Implementation of \code{\link[base]{ifelse}} for SpatRasters. This method allows for a concise approach to what can otherwise be achieved with a combination of \code{\link{classify}}, \code{\link{mask}}, and \code{\link{cover}}. 

\code{ifel} is an \code{R} equivalent to the \code{Con} method in ArcGIS (arcpy). 
}

\usage{
\S4method{ifel}{SpatRaster}(test, yes, no, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{test}{SpatRaster}
  \item{yes}{SpatRaster or numeric}
  \item{no}{SpatRaster or numeric}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  
}

\value{
SpatRaster
}

\examples{
r <- rast(nrows=5, ncols=5)
values(r) <- -10:14

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

z <- ifel(r > -2 & r < 2, 100, 0)

# nested expressions
y <- ifel(r > 1, 1, ifel(r < -1, -1, r))

k <- ifel(r > 0, r+10, ifel(r < 0, r-10, 3))
}

\keyword{spatial}

\name{sprc}

\docType{methods}

\alias{src}
\alias{src,ANY-method}
\alias{sprc}
\alias{sprc,missing-method}
\alias{sprc,list-method}
\alias{sprc,SpatRaster-method}


\title{Create a SpatRasterCollection}

\description{
Methods to create a SpatRasterCollection. This is an object to hold a collection (list) of SpatRaster objects. There are no restrictions on the similarity of the SpatRaster geometry. 

They can be used to combine several SpatRasters to be used with \code{\link{merge}} or \code{\link{mosaic}}
}

\usage{
\S4method{sprc}{SpatRaster}(x, ...) 

\S4method{sprc}{list}(x) 

\S4method{sprc}{missing}(x) 
}

\arguments{
\item{x}{SpatRaster, list with SpatRaster objects, or missing}
\item{...}{additional SpatRaster objects} 
}

\value{
SpatRasterCollection
}

\seealso{
\code{\link{sds}}
}

\examples{
x <- rast(xmin=-110, xmax=-50, ymin=40, ymax=70, ncols=60, nrows=30)
y <- rast(xmin=-80, xmax=-20, ymax=60, ymin=30)
res(y) <- res(x)
values(x) <- 1:ncell(x)
values(y) <- 1:ncell(y)

z <- sprc(x, y)
z
}


\keyword{methods}
\keyword{spatial}


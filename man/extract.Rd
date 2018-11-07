\name{extract}

\docType{methods}

\alias{extract}
\alias{extract,SpatRaster,SpatLayer-method}

\title{Extract values from a SpatRaster}

\description{
Extract values from a SpatRaster for the locations of a SpatVector
} 

\usage{
\S4method{extract}{SpatRaster,SpatLayer}(x, y, fun="", ...)
}

\arguments{
\item{x}{SpatRaster object}
\item{y}{SpatVector object}
\item{fun}{character. function to summarize the data. Currently ignored}
\item{...}{additional arguments. None implemented}
}

\value{data.frame}

\seealso{\link{values}}

\examples{
r <- rast(ncol=10, nrow=10)
values(r) <- 1:100
p <- vect(cbind(0,0))
#extract(r, p)
}

\keyword{methods}
\keyword{spatial}


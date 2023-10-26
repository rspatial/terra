\name{freq}

\docType{methods}

\alias{freq}
\alias{freq,SpatRaster-method}

\title{Frequency table}

\description{
Frequency table of the values of a SpatRaster. \code{NA}s are not counted unless \code{value=NA}. 

You can provide a SpatVector or additional SpatRaster to define zones for which to do tabulations.
}

\usage{
\S4method{freq}{SpatRaster}(x, digits=0, value=NULL, bylayer=TRUE, usenames=FALSE, zones=NULL, wide=FALSE)
}

\arguments{
  \item{x}{SpatRaster}
  \item{digits}{integer. Used for rounding the values before tabulation. Ignored if \code{NA}}
  \item{value}{numeric. An optional single value to only count the number of cells with that value. This value can be \code{NA}}
  \item{bylayer}{logical. If \code{TRUE} tabulation is done by layer}  
  \item{usenames}{logical. If \code{TRUE} layers are identified by their names instead of their numbers. Only relevant if \code{bylayer} is \code{TRUE}}
  \item{zones}{SpatRaster or SpatVector to define zones for which the tabulation should be done}
  \item{wide}{logical. Should the results by "wide" instead of "long"?}  
}

\value{
A \code{data.frame} with 3 columns (layer, value, count) unless \code{bylayer=FALSE} in which case a\code{data.frame} with two columns is returned (value, count). 
}


\examples{
r <- rast(nrows=10, ncols=10)
set.seed(2)
values(r) <- sample(5, ncell(r), replace=TRUE)

freq(r)

x <- c(r, r/3)
freq(x, bylayer=FALSE)
freq(x)

freq(x, digits=1)
freq(x, digits=-1)

freq(x, value=5)

}

\keyword{spatial}
\keyword{univar}

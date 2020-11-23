\name{sds}

\docType{methods}

\alias{sds}
\alias{sds,missing-method}
\alias{sds,character-method}
\alias{sds,list-method}
\alias{sds,SpatRaster-method}
\alias{[<-,SpatDataSet,numeric,missing-method}


\title{Create a SpatDataSet}

\description{
Methods to create a SpatDataSet. This is an object to hold "sub-datasets", each a SpatRaster that in most cases will have multiple layers. 

See \code{\link{describe_sds}} for getting information about sub-datasets in a file.
}

\usage{
\S4method{sds}{missing}(x, ...) 

\S4method{sds}{character}(x, ids=0, ...)

\S4method{sds}{SpatRaster}(x, ...) 

\S4method{sds}{list}(x, ...) 
}

\arguments{
\item{x}{character (filename), or SpatRaster, or list of SpatRaster objects, or missing. If multiple filenames are provided, it is attempted to make SpatRasters from these, and combine them into a SpatDataSet}
\item{ids}{optional. vector of integer subdataset ids. Ignored if the first value is not a postive integer}
\item{...}{additional arguments. Can be other \code{SpatRaster} objects if \code{x} is a \code{SpatRaster}} 
}


\value{
SpatDataSet
}

\seealso{
\code{\link{describe_sds}}
}

\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
x <- sds(s, s/2)
names(x) <- c("first", "second")
x
length(x)

# extract the second SpatRaster
x[2]
}


\keyword{methods}
\keyword{spatial}


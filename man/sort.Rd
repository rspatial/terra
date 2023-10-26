\name{sort}

\docType{methods}

\alias{sort}
\alias{sort,SpatRaster-method}
\alias{sort,SpatVector-method}
\alias{sort,data.frame-method}


\title{Sort a SpatRaster or SpatVector}

\description{
Sort the cell values of a SpatRaster across layers. You can also compute the sorting order. 

Or sort the records of SpatVector (or data.frame) by specifying the column number(s) or names(s) to sort on.
}

\usage{
\S4method{sort}{SpatRaster}(x, decreasing=FALSE, order=FALSE, filename="", ...)

\S4method{sort}{SpatVector}(x, v, decreasing=FALSE)
}

\arguments{
  \item{x}{SpatRaster}
  \item{decreasing}{logical. If \code{TRUE}, sorting is in decreasing order}
  \item{order}{logical. If \code{TRUE} the sorting order is returned instead of the sorted values}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
  \item{v}{character or numeric indicating the column(s) to sort on}
}


\value{
SpatRaster
}


\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
r <- c(r, r/2, r*2)
sort(r)

ord <- sort(r, order=TRUE)
# these two are the same 
ord[[1]]
which.min(r)
}

\keyword{spatial}


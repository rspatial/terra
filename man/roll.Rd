\name{roll}

\alias{roll}
\alias{roll,SpatRaster-method}
\alias{roll,numeric-method}

\title{Rolling (moving) functions}

\description{
Compute "rolling" or "moving" functions, such as the "rolling average"
}

\usage{
\S4method{roll}{SpatRaster}(x, n, fun=mean, type="around", circular=FALSE, na.rm=FALSE, filename="", ...) 

\S4method{roll}{numeric}(x, n, fun=mean, type="around", circular=FALSE, na.rm=FALSE) 
}
 

\arguments{
  \item{x}{SpatRaster or numeric}
  \item{n}{integer > 1. The size of the "window", that is, the number of sequential elements to use in the \code{fun}. This is normally an odd number}
  \item{fun}{a function like mean, min, max, sum}
  \item{type}{character. One of "around", "to", or "from". The choice indicates which values should be used in the computation. The focal element is always used. If \code{type} is "around", the other elements are before and after the focal element. Alternatively, you can select the elements preceding the focal element ("to") or those coming after it {"from"}. For example, when using \code{n=3} for element 5 of a vector; "around" used elements 4,5,6; "to" used elements 3,4,5, and "from" uses elements 5,6,7}
  \item{circular}{logical. If \code{TRUE}, the data are considered to have a circular nature (e.g. days or months of the year), such that there are no missing values before first or after the last value.}
  \item{na.rm}{logical. If \code{TRUE}, \code{NA} values should be ignored (by \code{fun})}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}


\value{
Same as \code{x}
}

\seealso{\code{\link{cumsum}}}

\examples{
## numeric
roll(1:12, 3, mean)
roll(1:12, 3, mean, "to")
roll(1:12, 3, mean, circular=TRUE)

## SpatVector
r <- rast(ncol=2, nrow=2, nlyr=10, vals=1)
r[1,2] = 2
r[2,2] = 4

roll(r, n=3, "sum", "from", na.rm=FALSE) |> values()
roll(r, n=3, "sum", "from", na.rm=TRUE) |> values()
roll(r, n=3, "sum", "from", circular=TRUE) |> values()

roll(r, n=3, "sum", "to", na.rm=TRUE) |> values()

roll(r, n=3, "sum", "around", circular=TRUE) |> values()
}

\keyword{spatial}

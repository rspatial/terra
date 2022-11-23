\name{roll}

\alias{roll}
\alias{roll,SpatRaster-method}
\alias{roll,numeric-method}

\title{Rolling (moving) functions}

\description{
Compute "rolling" or "moving" functions, such as the "rolling average"
}

\usage{
roll(x, n, fun=mean, type="around", circular=FALSE, na.rm=FALSE) 
}
 

\arguments{
  \item{x}{SpatRaster or numeric}
  \item{n}{Size of the 'window', i.e. the number of sequential elements to use in the function}
  \item{fun}{A function like mean, min, max, sum}
  \item{type}{Character. One of 'around', 'to', or 'from'. The choice indicates which values should be used in the computation. The focal element is always used. If \code{type} is 'around', the other elements are before and after the focal element. Alternatively, you can select the elements preceding the focal element ('to') or those coming after it {'from'}. For example, when using \code{n=3} for element 5 of a vector; 'around' used elements 4,5,6; 'to' used elements 3,4,5, and 'from' uses elements 5,6,7}
  \item{circular}{Logical. If \code{TRUE}, the data are considered to have a circular nature (e.g. months of the year), and the last elements in vector \code{x} are used in the computation of the moving function of the first element(s) of the vector, and the first elements are used in the computation of the moving function for the last element(s)}
  \item{na.rm}{Logical. If \code{TRUE}, \code{NA} values should be ingored (by \code{fun})}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}


\value{
Same as \code{x}
}

\seealso(\code{\link{cumsum}}

\examples{
## numeric
roll(1:12, 3, mean)
roll(1:12, 3, mean, 'to')
roll(1:12, 3, mean, 'from')
roll(1:12, 3, mean, circular=TRUE)

v <- c(0,1,2,3,3,3,3,4,4,4,5,5,6,7,7,8,9,NA)
roll(v, n=5)
roll(v, n=5, na.rm=TRUE)

## SpatVector
r <- rast(ncol=2, nrow=2, nlyr=10, vals=1)
r[1,2] = 2
r[2,2] = 4

roll(r, n=3, "sum", "from", na.rm=F) |> values()
roll(r, n=3, "sum", "from", na.rm=T) |> values()
roll(r, n=3, "sum", "from", circular=TRUE) |> values()

roll(r, n=3, "sum", "to", na.rm=T) |> values()

roll(r, n=3, "sum", "around", TRUE, na.rm=F) |> values()
}

\keyword{spatial}

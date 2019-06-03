\name{reclassify}

\docType{methods}

\alias{reclassify}
\alias{reclassify,SpatRaster-method}

\title{Reclassify}

\description{
Reclassify values of a SpatRaster. The function (re)classifies groups of values to other values. 

Reclassification can be based on ranges (from-to-becomes) or on specific values (is-becomes).

Reclassification is done with matrix \code{rcl}, in the row order of the reclassify table. Thus, if there are overlapping ranges or values, the first time a number is within a range determines the reclassification value.  
}

\usage{
\S4method{reclassify}{SpatRaster}(x, rcl, include.lowest=FALSE, right=TRUE, othersNA=FALSE, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{Raster* object}
  \item{rcl}{matrix for reclassification. This matrix must have 2 or 3 columns. If there are three columns, the first two columns are "from" "to" of the input values, and the third column "becomes" has the new value for that range. You can also supply a vector that is then coverted into a 3 column matrix with \code{byrow=TRUE}.  
  The two column matrix ("is", "becomes") can be useful for reclassifying integer values. In that case, the \code{right} argument is automatically set to \code{NA}}
 \item{include.lowest}{logical, indicating if a value equal to the lowest value in \code{rcl} (or highest value in the second column, for \code{right=FALSE}) should be included.}
 \item{right}{logical, indicating if the intervals should be closed on the right (and open on the left) or vice versa. The default is \code{TRUE}. A special case is to use right=NA. In this case both the left and right intervals are open}
 \item{othersNA}{logical. If \code{TRUE}, values that are not matched become \code{NA}. If \code{FALSE}, they retain their original value.} 
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}

 \item{...}{additional arguments. None implemented}  
}

\value{
SpatRaster
}

\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- runif(ncell(r)) 
# reclassify the values into three groups 
# all values >= 0 and <= 0.25 become 1, etc.
m <- c(0, 0.25, 1,  0.25, 0.5, 2,  0.5, 1, 3)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(r, rclmat)

# equivalent to
rc <- reclassify(r, c(-Inf,0.25,1, 0.25,0.5,2, 0.5,Inf,3))

# is becomes 
x <- round(r*5)
unique(x)
m <- rbind(c(1,100), c(2,200))
m
rcx1 <- reclassify(x, m)
unique(rcx1)
rcx2 <- reclassify(x, m, othersNA=TRUE)
unique(rcx2)

}

\keyword{spatial}

\name{classify}

\docType{methods}

\alias{classify}
\alias{classify,SpatRaster-method}

\title{classify}

\description{
Classify values of a SpatRaster. The function (re)classifies groups of values to other values. 

Classification can be based on ranges "from-to-becomes" or on specific values "is-becomes", or on "cuts".

With "from-to-becomes" or "is-becomes", classification is done with matrix \code{rcl}, in the row order of the matrix. Thus, if there are overlapping ranges or values, the first time a number is within a range determines the reclassification value. 

With "cuts" the values are sorted, so that the order in which they are provided does not matter. 
}

\usage{
\S4method{classify}{SpatRaster}(x, rcl, include.lowest=FALSE, right=TRUE, othersNA=FALSE, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{rcl}{matrix for classification. This matrix must have 1, 2 or 3 columns. If there are three columns, the first two columns are "from" "to" of the input values, and the third column "becomes" has the new value for that range. 
  The two column matrix ("is", "becomes") can be useful for re-classifying integer values. In that case, the \code{right} argument is automatically set to \code{NA}.
  A single column matrix or a vector is interpreted as a set of cuts. In that case the values are classified based on their location inbetween the cut-values}  
    
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

\note{
For model-based classification see \code{\link{predict}}
}

\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- (0:99)/99
# classify the values into three groups 
# all values >= 0 and <= 0.25 become 1, etc.
m <- c(0, 0.25, 1,  
       0.25, 0.5, 2,
	   0.5, 1, 3)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc1 <- classify(r, rclmat, include.lowest=TRUE)

# equivalent to
rc2 <- classify(r, c(0, 0.25, 0.5, 1), include.lowest=TRUE)

# is becomes 
x <- round(r*5)
unique(x)
m <- rbind(c(1,100), c(2,200))
m
rcx1 <- classify(x, m)
unique(rcx1)
rcx2 <- classify(x, m, othersNA=TRUE)
unique(rcx2)
}

\keyword{spatial}

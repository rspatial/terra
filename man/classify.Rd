\name{classify}

\docType{methods}

\alias{classify}
\alias{classify,SpatRaster-method}

\title{Classify (or reclassify) cell values}

\description{
Classify values of a SpatRaster. The function (re-)classifies groups of values to other values. 

The classification is done based on the argument \code{rcl}. You can classify ranges by specifying a three-column matrix "from-to-becomes" or change specific values by using a two-column matrix "is-becomes". You can also supply a vector with "cuts" or the "number of cuts".

With "from-to-becomes" or "is-becomes" classification is done in the row order of the matrix. Thus, if there are overlapping ranges or values, the first time a number is within a range determines the reclassification value. 

With "cuts" the values are sorted, so that the order in which they are provided does not matter. 
}

\usage{
\S4method{classify}{SpatRaster}(x, rcl, include.lowest=FALSE, right=TRUE, 
     others=NULL, brackets=TRUE, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{rcl}{matrix for classification. This matrix must have 1, 2 or 3 columns. If there are three columns, the first two columns are "from" "to" of the input values, and the third column "becomes" has the new value for that range. 
  
  The two column matrix ("is", "becomes") can be useful for classifying integer values. In that case, the arguments \code{right} and \code{include.lowest} are ignored.
  
  A single column matrix (or a vector) is interpreted as a set of cuts if there is more than one value. In that case the values are classified based on their location in-between the cut-values.  
  
  If a single number is provided, that is used to make that number of cuts, at equal intervals between the lowest and highest values of the SpatRaster.
  } 
  
  \item{include.lowest}{logical, indicating if a value equal to the lowest value in \code{rcl} (or highest value in the second column, for \code{right=FALSE}) should be included.}
  
 \item{right}{logical. If \code{TRUE}, the intervals are closed on the right (and open on the left). If \code{FALSE} they are open at the right and closed at the left. "open" means that the extreme value is *not* included in the interval. Thus, right-closed and left open is \code{(0,1] = {x | 0 < x <= 1}}. You can also close both sides with \code{right=NA}, that is only meaningful if you "from-to-becomes" classification with integers. For example to classify 1-5 -> 1, 6-10 -> 2, 11-15 -> 3. That may be easier to read and write than the equivalent 1-5 -> 1, 5-10 -> 2, 10-15 -> 3 with \code{right=TRUE} and \code{include.lowest=TRUE}}
 
 \item{others}{numeric. If not \code{NULL} all values that are not matched are set to this value. Otherwise they retain their original value.} 

 \item{brackets}{logical. If \code{TRUE}, intervals are have parenthesis or brackets around them to indicate whether they are open or closed. Only applies if \code{rcl} is a vector (or single column matrix)} 

 \item{filename}{character. Output filename}
 \item{...}{Additional arguments for writing files as in \code{\link{writeRaster}}}

}

\value{
SpatRaster
}

\seealso{\code{\link{subst}} for simpler from-to replacement}

\note{
classify works with the "raw" values of categorical rasters, ignoring the levels (labels, categories). To change the labels of categorical rasters, use \code{\link{subst}} instead.

For model-based classification see \code{\link{predict}}
}

\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- (0:99)/99

## from-to-becomes
# classify the values into three groups 
# all values >= 0 and <= 0.25 become 1, etc.
m <- c(0, 0.25, 1,
       0.25, 0.5, 2,
       0.5, 1, 3)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc1 <- classify(r, rclmat, include.lowest=TRUE)

## cuts
# equivalent to the above, but now a categorical SpatRaster is returned 
rc2 <- classify(r, c(0, 0.25, 0.5, 1), include.lowest=TRUE, brackets=TRUE)
freq(rc2)

## is-becomes 
x <- round(r*3)
unique(x)
# replace 0 with NA
y <- classify(x, cbind(0, NA))
unique(y)

# multiple replacements
m <- rbind(c(2, 200), c(3, 300))
m

rcx1 <- classify(x, m)
unique(rcx1)

rcx2 <- classify(x, m, others=NA)
unique(rcx2)
}

\keyword{spatial}


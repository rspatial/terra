\name{focal}

\alias{focal}
\alias{focal,SpatRaster-method}


\title{Focal values}

\description{
Calculate focal ("moving window") values for the neighborhood of focal cells using a matrix of weights, perhaps in combination with a function.
}

\usage{
\S4method{focal}{SpatRaster}(x, w=3, na.rm=TRUE, fillvalue=NA, fun="sum", filename="", format="", datatype="FLT4S", overwrite=FALSE, ...) 
}

\arguments{
\item{x}{SpatRaster}
  
\item{w}{matrix of weights (the moving window), e.g. a 3 by 3 matrix with values 1; see Details. The matrix does not need to be square, but the sides must be odd numbers. If you need even sides, you can add a column or row with weights of zero}

\item{na.rm}{logical. If \code{TRUE}, \code{NA} will be removed from focal computations. The result will only be \code{NA} if all focal cells are \code{NA}. Except for some special cases (weights of 1, functions like min, max, mean), using \code{na.rm=TRUE} is generally not a good idea in this function because it will unbalance the effect of the weights }

\item{fillvalue}{numeric. The value of the cells of the padded rows and columns}

  
\item{fun}{function (optional). The function fun should take multiple numbers, and return a single number. For example mean, modal, min or max. It should also accept a \code{na.rm} argument (or ignore it, e.g. as one of the 'dots' arguments. For example, \code{length} will fail, but \code{function(x, ...){na.omit(length(x))}} works. }

\item{filename}{Character. Output filename. Optional}
\item{format}{Character. Output file format. See \code{\link{writeRaster}}}
\item{datatype}{Character. Output file data type. See \code{\link{writeRaster}}} 
\item{overwrite}{Logical. If \code{TRUE}, \code{filename} is overwritten}
\item{...}{Additional arguments. None implemented}
}

\details{
\code{focal} uses a matrix of weights for the neighborhood of the focal cells. The default function is \code{sum}. It is computationally much more efficient to adjust the weights-matrix than to use another function through the \code{fun} argument. Thus while the following two statements are equivalent (if there are no \code{NA} values), the first one is faster than the second one:

\code{a <- focal(x, w=matrix(1/9, nc=3, nr=3))}

\code{b <- focal(x, w=matrix(1,3,3), fun=mean)}

There is, however, a difference if \code{NA} values are considered. One can use the \code{na.rm=TRUE} option which may make sense when using a function like \code{mean}. However, the results would be wrong when using a weights matrix.

Laplacian filter: \code{filter=matrix(c(0,1,0,1,-4,1,0,1,0), nrow=3)}

Sobel filter: \code{filter=matrix(c(1,2,1,0,0,0,-1,-2,-1) / 4, nrow=3)}
}

\value{SpatRaster object}


\examples{
r <- rast(ncols=36, nrows=18, xmn=0)
values(r) <- runif(ncell(r)) 

# 3x3 mean filter
# r <- focal(r, w=matrix(1/9,nrow=3,ncol=3), fun=mean) 
}



\keyword{spatial}


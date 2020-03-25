\name{focal}

\alias{focal}
\alias{focal,SpatRaster-method}


\title{Focal values}

\description{
Calculate focal ("moving window") values for the neighborhood of focal cells using a matrix of weights, perhaps in combination with a function.
}

\usage{
\S4method{focal}{SpatRaster}(x, w=3, na.rm=TRUE, fillvalue=NA, fun="sum", 
            filename="", overwrite=FALSE, wopt=list(), ...) 
}

\arguments{
\item{x}{SpatRaster}
  
\item{w}{window. The window can be defined as one (for a square) or two numbers (row, col); or with an odd-sized matrix with weights. See Details. If the window is defined as a matrix, \code{na.rm=FALSE} and \code{fun=sum}, and cannot be changed}

\item{na.rm}{logical. Should missing values be removed?}

\item{fillvalue}{numeric. The value of the cells of the padded rows and columns}

\item{fun}{function. The function fun should take multiple numbers, and return a single number. For example mean, modal, min or max. It should also accept a \code{na.rm} argument}


\item{filename}{character. Output filename. Optional}
\item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
\item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
\item{...}{additional arguments. None implemented}  
}

\details{
\code{focal} 

The window used must have odd dimensions. If you need even sides, you can make a matrix add a column or row with weights of zero. 

Example weight matrices

Laplacian filter: \code{filter=matrix(c(0,1,0,1,-4,1,0,1,0), nrow=3)}

Sobel filter: \code{filter=matrix(c(1,2,1,0,0,0,-1,-2,-1) / 4, nrow=3)}
}

\value{SpatRaster}


\examples{
r <- rast(ncols=10, nrows=10, ext(0, 10, 0, 10))
values(r) <- 1:ncell(r)

f <- focal(r, w=3, fun=function(x, ...)quantile(x, .25, ...), na.rm=TRUE) 

f <- focal(r, w=3, fun="mean") 

# the following two statements are equivalent:
a <- focal(r, w=matrix(1/9, nc=3, nr=3))
b <- focal(r, w=3, fun=mean, na.rm=FALSE)

# but this is different
d <- focal(r, w=3, fun=mean, na.rm=TRUE)
}


\keyword{spatial}


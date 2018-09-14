\name{focal}

\alias{focal}
\alias{focal,SpatRaster-method}


\title{Focal values}

\description{
Calculate focal ("moving window") values for the neighborhood of focal cells using a matrix of weights, perhaps in combination with a function.
}

\usage{
\S4method{focal}{SpatRaster}(x, w, fun, filename='', na.rm=FALSE, pad=FALSE, padValue=NA, NAonly=FALSE, ...) 
}

\arguments{
\item{x}{SpatRaster}
  
\item{w}{matrix of weights (the moving window), e.g. a 3 by 3 matrix with values 1; see Details. The matrix does not need to be square, but the sides must be odd numbers. If you need even sides, you can add a column or row with weights of zero}
  
\item{fun}{function (optional). The function fun should take multiple numbers, and return a single number. For example mean, modal, min or max. It should also accept a \code{na.rm} argument (or ignore it, e.g. as one of the 'dots' arguments. For example, \code{length} will fail, but \code{function(x, ...){na.omit(length(x))}} works. }

\item{filename}{character. Filename for a new raster (optional)}
 
\item{na.rm}{logical. If \code{TRUE}, \code{NA} will be removed from focal computations. The result will only be \code{NA} if all focal cells are \code{NA}. Except for some special cases (weights of 1, functions like min, max, mean), using \code{na.rm=TRUE} is generally not a good idea in this function because it will unbalance the effect of the weights }
 
\item{pad}{logical. If \code{TRUE}, additional 'virtual' rows and columns are padded to \code{x} such that there are no edge effects. This can be useful when a function needs to have access to the central cell of the filter}
  
\item{padValue}{numeric. The value of the cells of the padded rows and columns}
 
\item{NAonly}{logical. If \code{TRUE}, only cell values that are \code{NA} are replaced with the computed focal values}
 
\item{...}{Additional arguments as for \code{\link{writeRaster}}}
}

\details{
\code{focal} uses a matrix of weights for the neighborhood of the focal cells. The default function is \code{sum}. It is computationally much more efficient to adjust the weights-matrix than to use another function through the \code{fun} argument. Thus while the following two statements are equivalent (if there are no \code{NA} values), the first one is faster than the second one:

\code{a <- focal(x, w=matrix(1/9, nc=3, nr=3))}

\code{b <- focal(x, w=matrix(1,3,3), fun=mean)}

There is, however, a difference if \code{NA} values are considered. One can use the \code{na.rm=TRUE} option which may make sense when using a function like \code{mean}. However, the results would be wrong when using a weights matrix.

Laplacian filter: \code{filter=matrix(c(0,1,0,1,-4,1,0,1,0), nrow=3)}

Sobel filter: \code{filter=matrix(c(1,2,1,0,0,0,-1,-2,-1) / 4, nrow=3)}

see the \code{\link{focalWeight}} function to create distance based circular, rectangular, or Gaussian filters.
}

\value{
RasterLayer
}


\seealso{ \code{\link{focalWeight}} }

\examples{
r <- georst(ncols=36, nrows=18, xmn=0)
values(r) <- runif(ncell(r)) 

# 3x3 mean filter
r3 <- focal(r, w=matrix(1/9,nrow=3,ncol=3)) 

# 5x5 mean filter
r5 <- focal(r, w=matrix(1/25,nrow=5,ncol=5)) 

# Gaussian filter
gf <- focalWeight(r, 2, "Gauss")
rg <- focal(r, w=gf)

# The max value for the lower-rigth corner of a 3x3 matrix around a focal cell
f = matrix(c(0,0,0,0,1,1,0,1,1), nrow=3)
f
rm <- focal(r, w=f, fun=max)

# global lon/lat data: no 'edge effect' for the columns
xmin(r) <- -180
r3g <- focal(r, w=matrix(1/9,nrow=3,ncol=3)) 


\dontrun{
## focal can be used to create a cellular automaton

# Conway's Game of Life 
w <- matrix(c(1,1,1,1,0,1,1,1,1), nr=3,nc=3)
gameOfLife <- function(x) {
	f <- focal(x, w=w, pad=TRUE, padValue=0)
	# cells with less than two or more than three live neighbours die
	x[f<2 | f>3] <- 0
	# cells with three live neighbours become alive
	x[f==3] <- 1
	x
}

# simulation function
sim <- function(x, fun, n=100, pause=0.25) {
	for (i in 1:n) {
		x <- fun(x)
		plot(x, legend=FALSE, asp=NA, main=i)
		dev.flush()
		Sys.sleep(pause)
	}
	invisible(x)
}

# Gosper glider gun
m <- matrix(0, nc=48, nr=34)
m[c(40, 41, 74, 75, 380, 381, 382, 413, 417, 446, 452, 480, 
  486, 517, 549, 553, 584, 585, 586, 619, 718, 719, 720, 752, 
  753, 754, 785, 789, 852, 853, 857, 858, 1194, 1195, 1228, 1229)] <- 1
init <- georst(m)

# run the model
sim(init, gameOfLife, n=150, pause=0.05)
}
}


\keyword{spatial}


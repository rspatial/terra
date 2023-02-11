\name{focalCpp}

\alias{focalCpp}
\alias{focalCpp,SpatRaster-method}


\title{Compute focal values with an iterating C++ function}

\description{
Calculate focal values with a C++ function that iterates over cells to speed up computations by avoiding an R loop (with apply).

See \code{\link{focal}} for an easier to use method. 
}

\usage{
\S4method{focalCpp}{SpatRaster}(x, w=3, fun, ..., fillvalue=NA, 
		silent=TRUE, filename="", overwrite=FALSE, wopt=list())
}

\arguments{
\item{x}{SpatRaster}
\item{w}{window. The window can be defined as one (for a square) or two numbers (row, col); or with an odd-sized weights matrix. See the Details section in \code{\link{focal}}}

\item{fun}{\code{\link{cppFunction}} that iterates over cells. For C++ functions that operate on a single focal window, or for R functions use \code{\link{focal}} instead. The function must have at least three arguments. The first argument can have any name, but it must be a \code{Rcpp::NumericVector}, \code{Rcpp::IntegerVector} or a \code{std::vector<double>}. This is the container that receives the focal values. The other two arguments \code{ni} and \code{wi} must be of type \code{size_t}. \code{ni} represents the number of cells and \code{nw} represents the size of (number of elements in) the window}

\item{...}{additional arguments to \code{fun}}

\item{fillvalue}{numeric. The value of the cells in the virtual rows and columns outside of the raster}

\item{silent}{logical. If \code{TRUE} error messages are printed that may occur when trying \code{fun} to determine the length of the returned value. This can be useful in debugging a \code{fun} that does not work}

\item{filename}{character. Output filename}
\item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}

\item{wopt}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{SpatRaster}


\seealso{ \code{\link{focal}}, \code{\link{focalValues}} }


\examples{
\dontrun{
library(Rcpp)
cppFunction( 
	'NumericVector sum_and_multiply(NumericVector x, double m, size_t ni, size_t nw) {
		NumericVector out(ni);
		// loop over cells
		size_t start = 0;
		for (size_t i=0; i<ni; i++) {
			size_t end = start + nw;
			// compute something for a window
			double v = 0;
			// loop over the values of a window
			for (size_t j=start; j<end; j++) {
				v += x[j];
			}
			out[i] = v * m;
			start = end;
		}
		return out;
	}'
)

nr <- nc <- 10
r <- rast(ncols=nc, nrows=nr, ext= c(0, nc, 0, nr))
values(r) <- 1:ncell(r)

raw <- focalCpp(r, w=3, fun=sum_and_multiply, fillvalue=0, m=10)

# same as
f1 <- focal(r, w=3, fun=sum, fillvalue=0) *10
all(values(f1) == values(raw))

# and as
ffun <- function(x, m) { sum(x) * m }
f2 <- focal(r, w=3, fun=ffun, fillvalue=0, m=10)


# You can also use an R function with focalCpp but this
# is not recommended 

R_sm_iter <- function(x, m, ni, nw) {
	out <- NULL
	for (i in 1:ni) {
		start <- (i-1) * nw + 1
		out[i] <- sum(x[start:(start+nw-1)]) * m
	}
	out
}

fr <- focalCpp(r, w=3, fun=R_sm_iter, fillvalue=0, m=10)

}
}


\keyword{spatial}


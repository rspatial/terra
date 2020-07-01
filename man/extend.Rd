\name{extend}

\alias{extend}
\alias{extend,SpatRaster-method}
\alias{extend,SpatExtent-method}


\title{Extend}

\description{
Extend returns an SpatRaster object with a larger spatial extent. See \code{\link{crop}} if you (also) want to remove rows or columns. 

There is also an extend method for SpatExtent objects to enlarge a SpatExtent. You can also use algebraic notation to do that (see examples)
}

\usage{
\S4method{extend}{SpatRaster}(x, y, filename="", overwrite=FALSE, wopt=list(), ...) 

\S4method{extend}{SpatExtent}(x, y, ...) 
}

\arguments{
  \item{x}{SpatRaster or SpatExtent}
  \item{y}{If \code{x} is a SpatRaster, \code{y} should be a SpatExtent, or an object from which it can be extracted (such as SpatRaster and SpatVector objects). Alternatively, you can provide a numeric vector of length 2 indicating the number of rows and columns that need to be added (or a single number when the number of rows and columns is equal)
  
  If \code{x} is a SpatExtent, \code{y} should be a numeric vector of 1, 2, or 4 elements}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  

}

\value{
SpatRaster or SpatExtent
}

\seealso{\code{\link{crop}}, \code{\link{merge}}, \code{\link{ext}}}

\examples{
r <- rast(xmin=-150, xmax=-120, ymin=30, ymax=60, ncol=36, nrow=18)
values(r) <- 1:ncell(r)
e <- ext(-180, -100, 40, 70)
re <- extend(r, e)

# extend with a number of rows and columns (at each side)
re2 <- extend(r, c(2,10))

# SpatExtent object
e <- ext(r)
e
extend(e, 10)
extend(e, 10, -10, 0, 20)
}

\keyword{spatial}

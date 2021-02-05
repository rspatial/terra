\name{expand}

\alias{extend}
\alias{extend,SpatRaster-method}
\alias{expand}
\alias{expand,SpatRaster-method}
\alias{expand,SpatExtent-method}


\title{Expand}

\description{
Expand the extent of a SpatRaster. See \code{\link{crop}} if you (also) want to remove rows or columns. 

You can also enlarge a SpatExtent with this method, or with algebraic notation (see examples)
}

\usage{
\S4method{expand}{SpatRaster}(x, y, filename="", overwrite=FALSE, wopt=list(), ...) 

\S4method{expand}{SpatExtent}(x, y, ...) 
}

\arguments{
  \item{x}{SpatRaster or SpatExtent}
  \item{y}{If \code{x} is a SpatRaster, \code{y} should be a SpatExtent, or an object from which it can be extracted (such as SpatRaster and SpatVector objects). Alternatively, you can provide two postive integers indicating the number of rows and columns that need to be added at each side (or a single postive integer when the number of rows and columns is equal)
  
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
re <- expand(r, e)

# expand with a number of rows and columns (at each side)
re2 <- expand(r, c(2,10))

# SpatExtent
e <- ext(r)
e
expand(e, 10)
expand(e, 10, -10, 0, 20)
}

\keyword{spatial}

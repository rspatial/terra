\name{aggregate}

\docType{methods}

\alias{aggregate}
\alias{aggregate,SpatRaster-method}
\alias{aggregate,SpatVector-method}

\title{Aggregate raster cells}

\description{
Aggregate a SpatRaster to create a new SpatRaster with a lower resolution (larger cells). Aggregation groups rectangular areas to create larger cells. The value for the resulting cells is computed with a user-specified function. 

Or aggregate ("dissolve") a SpatVector. 
}

\usage{
\S4method{aggregate}{SpatRaster}(x, fact=2, fun="mean", ..., filename="", overwrite=FALSE, wopt=list())

\S4method{aggregate}{SpatVector}(x, by=NULL, sums=NULL, dissolve=TRUE, vars=NULL, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{fact}{positive integer. Aggregation factor expressed as number of cells in each direction (horizontally and vertically). Or two integers (horizontal and vertical aggregation factor) or three integers (when also aggregating over layers)}  
  \item{fun}{function used to aggregate values. Either an actual function, or a character specifying a C++ level implemented function: "mean", "max", "min", "median", "sum" and "modal"}
  \item{...}{additional arguments passed to \code{fun}, such as \code{na.rm=TRUE}}  
  \item{nodes}{positive integer. If \code{nodes > 1}, a 'parallel' package cluster with that many nodes is created. Ignored for C++ level implemented functions "mean", "max", "min", "median", "sum" and "modal"}  
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  
  \item{by}{character}
  \item{sums}{list}
  \item{dissolve}{logical}
  \item{vars}{character}
   
}


\details{
Aggregation starts at the upper-left end of a SpatRaster. If a division of the number of columns or rows with \code{factor} does not return an integer, the extent of the resulting SpatRaster will be somewhat larger then that of the original SpatRaster. For example, if an input SpatRaster has 100 columns, and \code{fact=12}, the output SpatRaster will have 9 columns and the maximum x coordinate of the output SpatRaster is also adjusted.
  
The function \code{fun} should take multiple numbers, and return a single number. For example \code{mean}, \code{modal}, \code{min} or \code{max}.

It should also accept a \code{na.rm} argument (or ignore it as one of the 'dots' arguments). 
}

\seealso{\code{\link{disaggregate}}}

\value{
SpatRaster
}


\examples{
r <- rast()
# aggregated SpatRaster, no values
ra <- aggregate(r, fact=10)

values(r) <- runif(ncell(r))
# aggregated raster, max of the values
ra <- aggregate(r, fact=10, fun=max)

# multiple layers
s <- c(r, r*2)
x <- aggregate(s, 2)
}

\keyword{methods}
\keyword{spatial}


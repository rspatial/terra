\name{aggregate}

\docType{methods}

\alias{aggregate}
\alias{aggregate,SpatRaster-method}
\alias{aggregate,SpatVector-method}

\title{Aggregate raster or vector data}

\description{
Aggregate a SpatRaster to create a new SpatRaster with a lower resolution (larger cells). Aggregation groups rectangular areas to create larger cells. The value for the resulting cells is computed with a user-specified function. 

You can also aggregate ("dissolve") a SpatVector. This either combines all geometries into one geometry, or it combines the geometries that have the same value for the variable(s) specified with argument \code{by}.
}

\usage{
\S4method{aggregate}{SpatRaster}(x, fact=2, fun="mean", ..., cores=1, filename="", overwrite=FALSE, wopt=list())

\S4method{aggregate}{SpatVector}(x, by=NULL, dissolve=TRUE, fun="mean", count=TRUE, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{fact}{positive integer. Aggregation factor expressed as number of cells in each direction (horizontally and vertically). Or two integers (horizontal and vertical aggregation factor) or three integers (when also aggregating over layers)}  
  \item{fun}{function used to aggregate values. Either an actual function, or for the following, their name: "mean", "max", "min", "median", "sum", "modal", "any", "all", "prod", "which.min", "which.max", "sd" (sample standard deviation) and "std" (population standard deviation)}
  \item{...}{additional arguments passed to \code{fun}, such as \code{na.rm=TRUE}}  
  \item{cores}{positive integer. If \code{cores > 1}, a 'parallel' package cluster with that many cores is created. Ignored for C++ level implemented functions that are listed under \code{fun}}  
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
  
  \item{by}{character. The variable(s) used to group the geometries}
  \item{dissolve}{logical. Should borders between aggregated geometries be dissolved?} 
  \item{count}{logical. If \code{TRUE} and \code{by} is not \code{NULL}, a variable "agg_n" is included that shows the number of input geometries for each output geometry}
}


\details{
Aggregation starts at the upper-left end of a SpatRaster. If a division of the number of columns or rows with \code{factor} does not return an integer, the extent of the resulting SpatRaster will be somewhat larger then that of the original SpatRaster. For example, if an input SpatRaster has 100 columns, and \code{fact=12}, the output SpatRaster will have 9 columns and the maximum x coordinate of the output SpatRaster is also adjusted.
  
The function \code{fun} should take multiple numbers, and return a single number. For example \code{mean}, \code{modal}, \code{min} or \code{max}.

It should also accept a \code{na.rm} argument (or ignore it as one of the 'dots' arguments). 
}

\seealso{\code{\link{disagg}} to disaggregate}

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
x <- aggregate(s, 20)


## SpatVector 
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
va <- aggregate(v, "ID_1")

plot(va, "NAME_1", lwd=5, plg=list(x="topright"), mar=rep(2,4))
lines(v, lwd=3, col="light gray")
lines(va)
text(v, "ID_1", halo=TRUE)
}

\keyword{methods}
\keyword{spatial}


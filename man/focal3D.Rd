\name{focal3D}

\alias{focal3D}
\alias{focal3D,SpatRaster-method}

\title{Three-dimensional focal values}

\description{
Calculate focal ("moving window") values for the three-dimensional neighborhood (window) of focal cells. See \code{\link{focal}} for two-dimensional focal computation.
}

\usage{
\S4method{focal3D}{SpatRaster}(x, w=3, fun=mean, ..., na.policy="all", fillvalue=NA, pad=FALSE, 
	padvalue=fillvalue, expand=FALSE, silent=TRUE, 
	filename="", overwrite=FALSE, wopt=list()) 
}

\arguments{
\item{x}{SpatRaster}
  
\item{w}{window. A rectangular prism (cuboid) defined by three numbers or by a three-dimensional array. The values are used as weights, and are usually zero, one, NA, or fractions. The window used must have odd dimensions. If you desire to use even sides, you can use an array, and pad the values with rows and/or columns that contain only \code{NA}s.
}

\item{fun}{function that takes multiple numbers, and returns one or multiple numbers for each focal area. For example mean, modal, min or max}

\item{...}{additional arguments passed to \code{fun} such as \code{na.rm}}

\item{na.policy}{character. Can be used to determine the cells of \code{x}, in the central layer, for which focal values should be computed. Must be one of "all" (compute for all cells), "only" (only for cells that are \code{NA}) or "omit" (skip cells that are \code{NA}). Note that the value of this argument does not affect which cells around each focal cell are included in the computations (use \code{na.rm=TRUE} to ignore cells that are \code{NA} in the computation of the focal value)}

\item{fillvalue}{numeric. The value of the cells in the virtual rows and columns outside of the raster}

\item{pad}{logical. Add virtual layers before the first and after the last layer}

\item{padvalue}{numeric. The value of the cells in the virtual layers}

\item{expand}{logical. Add virtual layers before the first or after the last layer that are the same as the first or last layers. If \code{TRUE}, arguments \code{pad} and \code{padvalue} are ignored}

\item{silent}{logical. If \code{TRUE} error messages are printed that may occur when trying \code{fun} to determine the length of the returned value. This can be useful in debugging a function passed to \code{fun} that does not work}

\item{filename}{character. Output filename}
\item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}

\item{wopt}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{SpatRaster}

\seealso{ \code{\link{focal}}}

\examples{
r <- rast(system.file("ex/logo.tif", package="terra"))   
x <- focal3D(r, c(5,5,3), na.rm=TRUE)

a <- array(c(0,1,0,1,1,1,0,1,0, rep(1,9), 0,1,0,1,1,1,0,1,0), c(3,3,3))
a[a==0] <- NA
z <- focal3D(r, a, na.rm=TRUE)
}

\keyword{spatial}


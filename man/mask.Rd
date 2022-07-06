\name{mask}

\docType{methods}

\alias{mask}
\alias{mask,SpatRaster,SpatRaster-method}
\alias{mask,SpatRaster,SpatVector-method}
\alias{mask,SpatRaster,sf-method}
\alias{mask,SpatVector,SpatVector-method}
\alias{mask,SpatVector,sf-method}

\title{Mask values in a SpatRaster or SpatVector}

\description{
If \code{x} is a \code{SpatRaster}: 
Create a new SpatRaster that has the same values as SpatRaster \code{x}, except for the cells that are \code{NA} (or another \code{maskvalue}) in another SpatRaster (the 'mask'), or not covered by a SpatVector. These cells become \code{NA} (or another \code{updatevalue}).

If \code{x} is a \code{SpatVector}: 
Select geometries of \code{x} that intersect, or not intersect, with the geometries of \code{y}.
} 

\usage{
\S4method{mask}{SpatRaster,SpatRaster}(x, mask, inverse=FALSE, maskvalues=NA, 
   updatevalue=NA, filename="", ...)
   
\S4method{mask}{SpatRaster,SpatVector}(x, mask, inverse=FALSE, updatevalue=NA,
	touches=TRUE, filename="", ...)

\S4method{mask}{SpatRaster,sf}(x, mask, inverse=FALSE, updatevalue=NA,
	touches=TRUE, filename="", ...)

\S4method{mask}{SpatVector,SpatVector}(x, mask, inverse=FALSE)
}


\arguments{
\item{x}{SpatRaster or SpatVector}
\item{mask}{SpatRaster or SpatVector or sf}
\item{inverse}{logical. If \code{TRUE}, areas on mask that are _not_ the \code{maskvalue} are masked}
\item{maskvalues}{numeric. The value(s) in \code{mask} that indicates the cells of \code{x} that should become \code{updatevalue} (default = \code{NA})}
\item{updatevalue}{numeric. The value that cells of \code{x} should become if they are not covered by \code{mask} (and not \code{NA})}
 \item{touches}{logical. If \code{TRUE}, all cells touched by lines or polygons will be masked, not just those on the line render path, or whose center point is within the polygon}
\item{filename}{character. Output filename}
\item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{SpatRaster}

\seealso{\link{crop}}


\examples{
r <- rast(ncols=10, nrows=10)
m <- rast(ncols=10, nrows=10)
values(r) <- 1:100
set.seed(1965)
x <- round(3 * runif(ncell(r)))
x[x==0] <- NA
values(m) <- x
mr <- mask(r, m)
}

\keyword{methods}
\keyword{spatial}


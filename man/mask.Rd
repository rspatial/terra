\name{mask}

\docType{methods}

\alias{mask}
\alias{mask,SpatRaster,SpatRaster-method}
\alias{mask,SpatRaster,SpatVector-method}

\title{Mask values in a SpatRaster}

\description{
Create a new SpatRaster that has the same values as SpatRaster \code{x}, except for the cells that are \code{NA} (or another \code{maskvalue}) in another SpatRaster (the 'mask'), or not covered by a SpatVector. These cells become \code{NA} (or another \code{updatevalue}).
} 

\usage{
\S4method{mask}{SpatRaster,SpatRaster}(x, mask, inverse=FALSE, maskvalues=NA, 
   updatevalue=NA, filename="", overwrite=FALSE, wopt=list(), ...)
   
\S4method{mask}{SpatRaster,SpatVector}(x, mask, inverse=FALSE, updatevalue=NA,
	touches=is.lines(mask), filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
\item{x}{SpatRaster}
\item{mask}{SpatRaster}
\item{inverse}{logical. If \code{TRUE}, areas on mask that are _not_ the \code{maskvalue} are masked}
\item{maskvalues}{numeric. The value(s) in \code{mask} that indicates the cells of \code{x} that should become \code{updatevalue} (default = \code{NA})}
\item{updatevalue}{numeric. The value that cells of \code{x} should become if they are not covered by \code{mask} (and not \code{NA})}
 \item{touches}{logical. If \code{TRUE}, all cells touched by lines or polygons will be masked, not just those on the line render path, or whose center point is within the polygon}
\item{filename}{character. Output filename. Optional}
\item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
\item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
\item{...}{additional arguments. None implemented}  

}

\value{SpatRaster}

\seealso{\link{crop}}


\examples{
r <- rast(ncol=10, nrow=10)
m <- rast(ncol=10, nrow=10)
values(r) <- 1:100
set.seed(1965)
x <- round(3 * runif(ncell(r)))
x[x==0] <- NA
values(m) <- x
mr <- mask(r, m)
}

\keyword{methods}
\keyword{spatial}


\name{window}

\docType{methods}

\alias{window}
\alias{window<-}
\alias{window,SpatRaster-method}
\alias{window<-,SpatRaster-method}

\title{Set a window}

\description{

Assign a window (area of interest) to a SpatRaster with a \code{SpatExtent}, or set it to \code{NULL} to remove the window. This is similar to \code{\link{crop}} without actually creating a new dataset. 

The window is intersect with the extent of the SpatRaster. It is envisioned that in a future version, the window may also go outside these boundaries.
}

\usage{
\S4method{window}{SpatRaster}(x)<-value

\S4method{window}{SpatRaster}(x)
}

\arguments{
  \item{x}{SpatRaster}
  \item{value}{SpatExtent}
}

\seealso{
 \code{\link{crop}}, \code{\link{extend}}
}


\value{
none for \code{window<-} and logical for \code{window}
}

\examples{ 
f <- system.file("ex/elev.tif", package="terra") 
r <- rast(f)
global(r, "mean", na.rm=TRUE)
e <- ext(c(5.9, 6,49.95, 50))

window(r) <- e
global(r, "mean", na.rm=TRUE)
r

x <- rast(f)
xe <- crop(x, e)
global(xe, "mean", na.rm=TRUE)

b <- c(xe, r)
window(b)
b

window(r) <- NULL
r
}

\keyword{methods}
\keyword{spatial}


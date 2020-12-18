\name{window}

\docType{methods}

\alias{window}
\alias{window<-}
\alias{window,SpatRaster-method}
\alias{window<-,SpatRaster-method}

\title{Set a window}

\description{

\bold{Experimental}: Assign a window (area of interest) to a SpatRaster with a \code{SpatExtent}, or set it to \code{NULL} to remove the window. This is similar to \code{\link{crop}} and/or \code{\link{expand}}, without actually creating a new dataset.
}

\usage{
\S4method{window}{SpatRaster}(x, ...)<-value

\S4method{window}{SpatRaster}(x, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{value}{SpatExtent}
  \item{...}{additional arguments. None implemented}  
}

\seealso{
 \code{\link{crop}}, \code{\link{extend}}
}


\value{
none for \code{window<-} and logical for \code{window}
}

\examples{ 
f <- system.file("ex/test.tif", package="terra") 
r <- rast(f)
global(r, "mean", na.rm=TRUE)
e <- ext(c(179680, 180176, 331905, 332304))

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


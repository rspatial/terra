\name{Summary-methods}

\docType{methods}

\alias{Summary-methods}
\alias{mean,SpatRaster-method}

\title{ Summary methods }

\description{

Summary method compute cell (pixel) level statistics across layers, and most return a single layer (the exception is \code{range}). 

The following summary methods are available for SpatRaster objects: \code{mean, max, min, range, prod, sum, any, all}

See \code{\link[terra]{reduce}} to compute statistics that are not included here (e.g. median) or any other custom functions.
 
Because generic functions are used, the method applied is chosen based on the first argument: '\code{x}'. This means that if \code{r} is a SpatRaster, \code{mean(r, 5)} will work, but \code{mean(5, r)} will not work.
  
To summarize the values of all cells by layer see \code{\link[terra]{layerStats}}
}

\value{SpatRaster}

\examples{
set.seed(0)
r <- rast(nrow=10, ncol=10, nlyr=3)
values(r) <- runif(length(r))

x <- max(r)
x <- mean(r)
#r <- range(r, 1.2)

#s <- sum(r, r[[2]], 3)

}

\keyword{methods}
\keyword{spatial}


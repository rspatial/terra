\name{local}

\docType{methods}

\alias{Summary-methods}

\alias{Summary,SpatRaster-method}
\alias{mean,SpatRaster-method}

\title{Local statistics}

\description{

Compute cell (pixel) level "local" statistics across layers, and most return a single layer (the exception is \code{range}). 

The following summary methods are available for SpatRaster: \code{mean, max, min, range, prod, sum, any, all}

See \code{reduce} to compute statistics that are not included here (e.g. median) or any other custom functions.
 
Because generic functions are used, the method applied is chosen based on the first argument: '\code{x}'. This means that if \code{r} is a SpatRaster, \code{mean(r, 5)} will work, but \code{mean(5, r)} will not work.
  
To summarize the values of all cells by layer see \code{layerStats}.
}


\value{SpatRaster}

\seealso{ \code{\link{Math-methods}} }


\examples{
set.seed(0)
r <- rast(nrow=10, ncol=10, nlyr=3)
values(r) <- runif(size(r))

x <- mean(r)
x <- sum(r, r[[2]], 3)
max(r)
max(r, 0.5)
}

\keyword{methods}
\keyword{spatial}


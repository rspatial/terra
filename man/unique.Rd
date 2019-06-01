\name{unique}

\docType{methods}

\alias{unique}
\alias{unique,SpatRaster,ANY-method}

\title{Unique values}

\description{
This function returns the unique values in a SpatRaster.
}

\usage{
\S4method{unique}{SpatRaster,missing}(x, incomparables=0, ...) 
}

\arguments{
  \item{x}{SpatRaster}
  \item{incomparables}{the number of decimals to use in rounding the values before evaluating their uniquecess}
  \item{...}{additional arguments. none implemented}  
}


\value{
vector or matrix
}



\examples{
r <- rast(ncol=5, nrow=5)
values(r) <- c(1:10, 1:10, 1:5)
unique(r)
s <- c(r, r/2)
u <- unique(s)
#do.call(rbind, u)
#lapply(u, unique)
}

\keyword{spatial}

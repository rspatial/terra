\name{catalyze}

\docType{methods}

\alias{as.numeric}
\alias{as.numeric,SpatRaster-method}
\alias{catalyze}
\alias{catalyze,SpatRaster-method}
  
  
\title{Factors to numeric}

\description{

Change a categorical layer into one or more numerical layers. With \code{as.numeric} you can transfer the active category values to cell values in a non-categorical SpatRaster. \code{catalyze} creates new layers for each category.
}

\usage{
\S4method{as.numeric}{SpatRaster}(x, index=NULL, filename="", ...)

\S4method{catalyze}{SpatRaster}(x, filename="", ...)
}


\arguments{
  \item{x}{SpatRaster}
  \item{index}{positive integer or category indicating the category to use. If \code{NULL} the active category is used}
  \item{filename}{character. Output filename}  
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}  
}  


\value{
SpatRaster
}

\seealso{\code{\link{activeCat}}, \code{\link{cats}}}


\examples{
set.seed(0)
r <- rast(nrows=10, ncols=10)
values(r) <- sample(3, ncell(r), replace=TRUE) + 10
d <- data.frame(id=11:13, cover=c("forest", "water", "urban"), letters=letters[1:3], value=10:12)
levels(r) <- d
catalyze(r)

activeCat(r) <- 3
as.numeric(r)
}

\keyword{methods}
\keyword{spatial}


\name{activeCat}

\docType{methods}

\alias{activeCat}
\alias{activeCat,SpatRaster-method}
\alias{activeCat<-}
\alias{activeCat<-,SpatRaster-method}
  
  
\title{Active category}

\description{
Get or set the active category of a multi-categorical SpatRaster layer
}

\usage{
\S4method{activeCat}{SpatRaster}(x, layer=1)
\S4method{activeCat}{SpatRaster}(x, layer=1)<-value
}


\arguments{
  \item{x}{SpatRaster}
  \item{layer}{positive integer, the layer number or name}
  \item{value}{a data.frame (ID, category)}
}  


\value{
integer
}

\seealso{\code{\link{catalyze}}, \code{\link{cats}}}


\examples{
set.seed(0)
r <- rast(nrows=10, ncols=10)
values(r) <- sample(3, ncell(r), replace=TRUE) + 10
d <- data.frame(id=11:13, cover=c("forest", "water", "urban"), letters=letters[1:3], value=10:12)
levels(r) <- d

activeCat(r)
activeCat(r) <- 3
activeCat(r)
}

\keyword{methods}
\keyword{spatial}


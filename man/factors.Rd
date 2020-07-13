\name{factors}

\docType{methods}

\alias{is.factor}
\alias{is.factor,SpatRaster-method}
\alias{as.factor}
\alias{as.factor,SpatRaster-method}
\alias{levels}
\alias{levels,SpatRaster-method}
\alias{levels<-}
\alias{levels<-,SpatRaster-method}
\alias{rats}
\alias{rats,SpatRaster-method}
\alias{setRat}
  
\title{Factors}

\description{

These functions allow for defining a SpatRaster layer as a categorical variable. The cell values are an index (id), whereas the actual values are stored seperately, in a table (sometimes called "Raster Attribute Table"). This table is a data.frame. The first column in the RAT ("ID") has the unique cell values of the layer; this column should normally not be changed. The other columns can be of any basic type (factor, character, integer, numeric or logical). 

Function 'levels' returns the RAT for inspection. It can be modified and set using \code{levels(x) <- value}.

\code{as.factor} and \code{ratify} create a layer with a RAT table. \code{deratify} creates a single layer for a (or each) variable in the RAT table. 
}

\usage{
is.factor(x)
as.factor(x)
levels(x)
\S4method{rats}{SpatRaster}(x)
setRat(x, layer, rat)
}


\arguments{
  \item{x}{SpatRaster}
  \item{layer}{layer number}
  \item{rat}{data.frame}
}  


\value{
SpatRaster; list (levels); boolean (is.factor)
}

\examples{
set.seed(0)
r <- rast(nrow=10, ncol=10)
values(r) <- sample(3, ncell(r), replace=TRUE) + 2
is.factor(r)
cls <- c("forest", "water", "urban")
levels(r) <- cls



f <- as.factor(r)
is.factor(f)
}

\keyword{methods}
\keyword{spatial}


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
\alias{rats<-}
\alias{rats<-,SpatRaster-method}
\alias{cats}
\alias{cats,SpatRaster-method}
\alias{cats<-}
\alias{cats<-,SpatRaster-method}
  
\title{Factors}

\description{

These functions allow for defining a SpatRaster layer as a categorical variable (factor). These can be inspected and set with \code{levels} and \code{cats}. These methods are the same except that with \code{levels<-} you can only set the categories of the first layer. Like factors in R, categories are stored as indices (integers) that have an associated label. In a SpatRaster, the index starts at 0, and cannot be larger than 255. \code{as.factor} create a layer with a RAT table. \code{deratify} creates a single layer for a (or each) variable in the RAT table. 


The "Raster Attribute Table" (RAT) is similar in concept. A main difference with categries is that there can be many labels (variables) assocuated with one index. The RAT is a data.frame; the first column ("ID") has the unique cell values of the layer.
}

\usage{
\S4method{is.factor}{SpatRaster}(x)

\S4method{as.factor}{SpatRaster}(x)

\S4method{levels}{SpatRaster}(x)

\S4method{levels}{SpatRaster}(x)<-value

\S4method{cats}{SpatRaster}(x)

\S4method{cats}{SpatRaster}(x, layer=1)<-value

\S4method{rats}{SpatRaster}(x)

\S4method{rats}{SpatRaster}(x, layer=1)<-value
}


\arguments{
  \item{x}{SpatRaster}
  \item{layer}{positive integer, the layer number or name}
  \item{value}{data.frame}
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


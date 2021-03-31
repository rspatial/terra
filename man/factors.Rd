\name{factors}

\docType{methods}

\alias{is.factor}
\alias{is.factor,SpatRaster-method}
\alias{levels}
\alias{levels,SpatRaster-method}
\alias{levels<-}
\alias{levels<-,SpatRaster-method}
\alias{cats}
\alias{cats,SpatRaster-method}
\alias{setCats}
\alias{setCats,SpatRaster-method}
  
\title{Factors}

\description{

A SpatRaster layer can be a categorical variable (factor). These can be inspected and set with \code{levels} and \code{cats}. With \code{levels<-} you can set the categories of the first layer and with \code{setCats} you can set it for any layer. Like factors in R, categories are stored as indices (integers) that have an associated label. In a SpatRaster, the index starts at 0, and cannot be larger than 255. 

There may be multiple categories (attributes) associated with a single layer. These can be set as a data.frame with \code{setCats} and they are returned with \code{cats}.
}

\usage{
\S4method{is.factor}{SpatRaster}(x)

\S4method{levels}{SpatRaster}(x)

\S4method{levels}{SpatRaster}(x)<-value

\S4method{cats}{SpatRaster}(x)

\S4method{setCats}{SpatRaster}(x, layer=1, value, index=1)

}


\arguments{
  \item{x}{SpatRaster}
  \item{layer}{positive integer, the layer number or name}
  \item{value}{ a characte vectors with new labels or a data.frame}
  \item{index}{positive integer, indicating the column in \code{value} to be used as the category (factor labels)}
}  


\value{
SpatRaster; list (e.g., levels); boolean (is.factor)
}

\examples{
set.seed(0)
r <- rast(nrow=10, ncol=10)
values(r) <- sample(3, ncell(r), replace=TRUE) - 1
is.factor(r)

cls <- c("forest", "water", "urban")
levels(r) <- cls
names(r) <- "land cover"
is.factor(r)
r

plot(r, col=c("green", "blue", "light gray"))
text(r, digits=3, cex=.75, halo=TRUE)

#p <- as.polygons(r)
#plot(p, "land cover", col=c("green", "blue", "light gray"))
}

\keyword{methods}
\keyword{spatial}


\name{names}

\alias{name}
\alias{name<-}
\alias{names}
\alias{names<-}

\alias{names,SpatRaster-method}
\alias{names<-,SpatRaster-method}

\alias{names,SpatRasterDataset-method}
\alias{names<-,SpatRasterDataset-method}

\alias{names,SpatRasterCollection-method}
\alias{names<-,SpatRasterCollection-method}

\alias{names,SpatVector-method}
\alias{names,SpatVectorCollection-method}
\alias{names,SpatVectorProxy-method}
\alias{names<-,SpatVector-method}
\alias{names<-,SpatVectorCollection-method}
  
\title{Names of Spat* objects}

\description{
Get or set the names of the layers of a SpatRaster or the attributes of a SpatVector. 

See \code{\link{set.names}} for in-place setting of names.
}

\usage{
\S4method{names}{SpatRaster}(x)

\S4method{names}{SpatRaster}(x)<-value

\S4method{names}{SpatRasterDataset}(x)

\S4method{names}{SpatRasterDataset}(x)<-value

\S4method{names}{SpatVector}(x)

\S4method{names}{SpatVector}(x)<-value
}

\arguments{
  \item{x}{SpatRaster, SpatRasterDataset, or SpatVector}
  \item{value}{character (vector)}
}

\value{
character 
}

\note{
terra enforces neither unique nor valid names. See \code{\link{make.unique}} to create unique names and \code{\link{make.names}} to make syntactically valid names.
}

\examples{
s <- rast(ncols=5, nrows=5, nlyrs=3)
nlyr(s)
names(s)
names(s) <- c("a", "b", "c")
names(s)

# SpatVector names
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
names(v)
names(v) <- paste0(substr(names(v), 1, 2), "_", 1:ncol(v))
names(v)
}

\keyword{spatial}


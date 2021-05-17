\name{names}

\alias{name}
\alias{name<-}
\alias{names}
\alias{names<-}
\alias{varnames}
\alias{varnames<-}
\alias{longnames}
\alias{longnames<-}

\alias{names,SpatRaster-method}
\alias{names<-,SpatRaster-method}
\alias{varnames,SpatRaster-method}
\alias{varnames<-,SpatRaster-method}

\alias{longnames,SpatRaster-method}
\alias{longnames<-,SpatRaster-method}

\alias{names,SpatRasterDataset-method}
\alias{names<-,SpatRasterDataset-method}
\alias{varnames,SpatRasterDataset-method}
\alias{varnames<-,SpatRasterDataset-method}
\alias{longnames,SpatRasterDataset-method}
\alias{longnames<-,SpatRasterDataset-method}

\alias{names,SpatVector-method}
\alias{names<-,SpatVector-method}
  
\title{Names of Spat* objects}

\description{
Get or set the names of the layers of a SpatRaster or the attributes of a SpatVector. With \code{lnames} you can get or set the "long names" of a SpatRaster or SpatRasterDataset. 

For a SpatRaster, you can also get/set a variable name or long name (one per data source)
}

\usage{
\S4method{names}{SpatRaster}(x)

\S4method{names}{SpatRaster}(x)<-value

\S4method{varnames}{SpatRaster}(x)

\S4method{varnames}{SpatRaster}(x)<-value

\S4method{longnames}{SpatRaster}(x)

\S4method{longnames}{SpatRaster}(x)<-value


\S4method{names}{SpatRasterDataset}(x)

\S4method{names}{SpatRasterDataset}(x)<-value

\S4method{varnames}{SpatRasterDataset}(x)

\S4method{varnames}{SpatRasterDataset}(x)<-value

\S4method{longnames}{SpatRasterDataset}(x)

\S4method{longnames}{SpatRasterDataset}(x)<-value

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
terra enforces neither unique nor valid names. See \code{\link{make.unique}} to create unique names and \code{{make.names}} to make syntactically valid names.
}

\examples{
s <- rast(ncols=5, nrows=5, nlyrs=3)
nlyr(s)
names(s)
names(s) <- c("a", "b", "c")
names(s)

# space is not valid
names(s)[2] <- "hello world"
names(s)

# two invalid names
names(s) <- c("a", "  a  ", "3")
names(s)

# SpatVector names
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
names(v)
names(v) <- paste0(substr(names(v), 1, 2), "_", 1:ncol(v))
names(v)
}

\keyword{spatial}


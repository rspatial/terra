\name{varnames}

\alias{varnames}
\alias{varnames<-}
\alias{longnames}
\alias{longnames<-}

\alias{varnames,SpatRaster-method}
\alias{varnames<-,SpatRaster-method}

\alias{longnames,SpatRaster-method}
\alias{longnames<-,SpatRaster-method}

\alias{varnames,SpatRasterDataset-method}
\alias{varnames<-,SpatRasterDataset-method}
\alias{longnames,SpatRasterDataset-method}
\alias{longnames<-,SpatRasterDataset-method}

 
\title{variable and long variable names}

\description{
Set or get names for each dataset (variable) in a SpatRasterDataset. 

Each SpatRaster _data source_ can also have a variable name and a long variable name. They are set when reading a file with possibly multiple sub-datasets (e.g. netcdf or hdf5 format) into a single SpatRaster. Each sub-datset is a seperate "data-source" in the SpatRaster. Note that newly created or derived SpatRasters always have a single variable (data source), and therefore the variable names are lost when processing a multi-variable SpatRaster. Thus the variable names are mostly useful to understand a SpatRaster created from some files and for managing SpatRasterDatasets.

See \code{link{names}} for the more commonly used _layer_ names.

}

\usage{
\S4method{varnames}{SpatRaster}(x)

\S4method{varnames}{SpatRaster}(x)<-value

\S4method{longnames}{SpatRaster}(x)

\S4method{longnames}{SpatRaster}(x)<-value

\S4method{varnames}{SpatRasterDataset}(x)

\S4method{varnames}{SpatRasterDataset}(x)<-value

\S4method{longnames}{SpatRasterDataset}(x)

\S4method{longnames}{SpatRasterDataset}(x)<-value
}

\arguments{
  \item{x}{SpatRaster, SpatRasterDataset}
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
names(s) <- c("a", "b", "c")
x <- sds(s, s)
varnames(x) <- c("one", "two")
x
}

\keyword{spatial}


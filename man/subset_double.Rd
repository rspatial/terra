\name{subset_double}


\alias{[[}

\alias{[[,SpatRaster,numeric,missing-method}
\alias{[[,SpatRaster,logical,missing-method}
\alias{[[,SpatRaster,character,missing-method}
\alias{[[,SpatRaster,ANY,missing-method}

\alias{[[,SpatRasterDataset,ANY,ANY-method}

\alias{[[,SpatVector,numeric,missing-method}
\alias{[[,SpatVector,logical,missing-method}
\alias{[[,SpatVector,character,missing-method}



\title{Subset a SpatRaster or a SpatVector}

\description{
Select a subset of layers from a SpatRaster or select a subset of records (row) and/or variables (columns) from a SpatVector.  
}

\usage{
\S4method{[[}{SpatRaster,numeric,missing}(x, i, j)

\S4method{[[}{SpatRasterDataset,ANY,ANY}(x, i, j, drop=TRUE)

\S4method{[[}{SpatVector,numeric,missing}(x, i, j, drop=FALSE)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{i}{
  if \code{x} is a \code{SpatRaster}: integer, logical, or character to select layers

  if \code{x} is a \code{SpatVector}: integer, logical, or character to select variables
  }
  \item{j}{missing, or, for SpatRasterDataset only, numeric}

  \item{drop}{logical. If \code{TRUE}, the geometries will be dropped, and a data.frame is returned}  
}

\value{
if \code{x} is a \code{SpatRaster} or \code{SpatRasterDataset}: SpatRaster

if \code{x} is a \code{SpatVector}: a \code{data.frame}.
}


\seealso{\code{\link{subset}}, \code{\link{$}}, \code{\link{[}}, \code{\link{extract}}}

\examples{
### SpatRaster
s <- rast(system.file("ex/logo.tif", package="terra"))   
s[[ 1:2 ]]

s[[c("red", "green")]]

# expression based (partial) matching of names with single brackets
s["re"]
s["^re"]

# does not with double brackets
# s[["re"]]


### SpatVector

v <- vect(system.file("ex/lux.shp", package="terra"))
v[[2:3]]

# to keep the geometry use
v[,2:3]
}

\keyword{ spatial }

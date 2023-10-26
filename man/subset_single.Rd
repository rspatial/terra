\name{subset_single}

\alias{[}

\alias{[,SpatExtent,missing,missing-method}
\alias{[,SpatExtent,numeric,missing-method}

\alias{[,SpatRaster,ANY,ANY,ANY-method}
\alias{[,SpatRaster,ANY,ANY-method}

\alias{[,SpatRasterDataset,logical,missing-method}
\alias{[,SpatRasterDataset,character,missing-method}

\alias{[,SpatRasterDataset,numeric,missing-method}
\alias{[,SpatRasterDataset,numeric,numeric-method}
\alias{[,SpatRasterDataset,numeric,logical-method}

\alias{[,SpatRasterDataset,missing,numeric-method}
\alias{[,SpatRasterDataset,missing,logical-method}

\alias{[,SpatRasterCollection,numeric,missing-method}

\alias{[,SpatVector,missing,missing-method}
\alias{[,SpatVector,missing,numeric-method}
\alias{[,SpatVector,missing,character-method}
\alias{[,SpatVector,missing,logical-method}
\alias{[,SpatVector,character,missing-method}
\alias{[,SpatVector,numeric,missing-method}
\alias{[,SpatVector,numeric,numeric-method}
\alias{[,SpatVector,numeric,character-method}
\alias{[,SpatVector,numeric,logical-method}
\alias{[,SpatVector,logical,missing-method}
\alias{[,SpatVector,logical,character-method}
\alias{[,SpatVector,logical,numeric-method}
\alias{[,SpatVector,logical,logical-method}
\alias{[,SpatVector,data.frame,missing-method}
\alias{[,SpatVector,data.frame,ANY-method}
\alias{[,SpatVector,matrix,missing-method}

\alias{[,SpatVector,SpatVector,missing-method}
\alias{[,SpatVector,SpatExtent,missing-method}
\alias{[,SpatRaster,SpatVector,missing-method}
\alias{[,SpatRaster,missing,missing-method}
\alias{[,SpatRaster,numeric,missing-method}
\alias{[,SpatRaster,missing,numeric-method}
\alias{[,SpatRaster,numeric,numeric-method}
\alias{[,SpatRaster,data.frame,missing-method}
\alias{[,SpatRaster,matrix,missing-method}
\alias{[,SpatRaster,SpatRaster,missing-method}
\alias{[,SpatRaster,SpatExtent,missing-method}


\title{Extract values from a SpatRaster, SpatVector or SpatExtent}

\description{
Extract values from a SpatRaster; a subset of records (row) and/or variables (columns) from a SpatVector; or a number from a SpatExtent.  

You can use indices (row, column, layer or cell numbers) to extract. You can also use other Spat* objects. 
}

\usage{
\S4method{[}{SpatRaster,ANY,ANY,ANY}(x, i, j, k)

\S4method{[}{SpatVector,numeric,numeric}(x, i, j, drop=FALSE)

\S4method{[}{SpatVector,SpatVector,missing}(x, i, j)

\S4method{[}{SpatExtent,numeric,missing}(x, i, j)

}

\arguments{
  \item{x}{SpatRaster, SpatVector or SpatExtent}
  \item{i}{
  if \code{x} is a \code{SpatRaster}: numeric, logical or missing to select rows or, if \code{j} is missing, to select cells numbers.
  
  if \code{x} is a \code{SpatVector}: numeric or missing to select rows. 
  if \code{i} is another SpatVector: get a new SpatVector with the geometries that intersect.

  if \code{x} is a \code{SpatExtent}: integer between 1 and 4.
  }
  \item{j}{numeric, logical, or missing to select columns}

  \item{k}{numeric, character, or missing to select layers}

  \item{drop}{logical. If \code{FALSE} an object of the same class as \code{x} is returned}
}

\value{
numeric if \code{x} is a SpatExtent. 
Same as \code{x} if \code{drop=FALSE}. Otherwise a data.frame
}

\seealso{\code{\link{extract}}, \code{\link{subset}}, \code{\link{$}}, \code{\link{[[}}}

\examples{
### SpatRaster
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
r[3638]
rowColFromCell(r, 2638)
r[39, 28]
x <- r[39:40, 28:29, drop=FALSE] 
as.matrix(x, wide=TRUE)

### SpatVector

v <- vect(system.file("ex/lux.shp", package="terra"))
v[2:3,]
v[1:2, 2:3]
v[1:2, 2:3, drop=TRUE]
}

\keyword{ spatial }

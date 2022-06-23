\name{subset}

\alias{subset}
\alias{subset,SpatRaster-method}
\alias{[[,SpatRaster,numeric,missing-method}
\alias{[[,SpatRaster,logical,missing-method}
\alias{[[,SpatRaster,character,missing-method}
\alias{[,SpatRaster,character,missing-method}
\alias{[,SpatRasterDataset,numeric,missing-method}
\alias{[,SpatRasterDataset,numeric,numeric-method}
\alias{[,SpatRasterDataset,logical,missing-method}
\alias{[,SpatRasterDataset,character,missing-method}
\alias{[[,SpatRasterDataset,ANY,ANY-method}

\alias{[,SpatRasterCollection,numeric,missing-method}

\alias{$,SpatRaster-method}
\alias{$,SpatRasterDataset-method}

\alias{subset,SpatVector-method}
\alias{[,SpatVector,missing,missing-method}
\alias{[,SpatVector,missing,numeric-method}
\alias{[,SpatVector,missing,character-method}
\alias{[,SpatVector,missing,logical-method}
\alias{[,SpatVector,numeric,missing-method}
\alias{[,SpatVector,numeric,numeric-method}
\alias{[,SpatVector,numeric,character-method}
\alias{[,SpatVector,numeric,logical-method}
\alias{[,SpatVector,logical,missing-method}
\alias{[,SpatVector,logical,character-method}
\alias{[,SpatVector,logical,numeric-method}
\alias{[,SpatVector,logical,logical-method}
\alias{[,SpatVector,data.frame,missing-method}
\alias{[,SpatVector,matrix,missing-method}



\title{Subset a SpatRaster or a SpatVector}

\description{
Select a subset of layers from a SpatRaster or select a subset of records (row) and/or variables (columns) from a SpatVector.  
}

\usage{
\S4method{subset}{SpatRaster}(x, subset, negate=FALSE, NSE=FALSE, filename="", overwrite=FALSE, ...)

\S4method{subset}{SpatVector}(x, subset, select, drop=FALSE, NSE=TRUE)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{subset}{
  if \code{x} is a \code{SpatRaster}: integer or character to select layers

  if \code{x} is a \code{SpatVector}: logical expression indicating the rows to keep (missing values are taken as FALSE)
  }
  \item{select}{expression, indicating columns to select} 
  \item{negate}{logical. If \code{TRUE} all layers that are \bold{not} in the subset are selected}
  \item{NSE}{logical. If \code{TRUE}, non-standard evaluation (the use of unquoted variable names) is allowed. Set this to \code{FALSE} when calling \code{subset} from a function} 
  \item{drop}{logical. If \code{TRUE}, the geometries will be dropped, and a data.frame is returned}  
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
if \code{x} is a \code{SpatRaster}: SpatRaster

if \code{x} is a \code{SpatVector}: SpatVector or, if \code{drop=TRUE}, a \code{data.frame}.
}


\examples{
### SpatRaster
s <- rast(system.file("ex/logo.tif", package="terra"))   
subset(s, 2:3)
subset(s, c(3,2,3,1))
#equivalent to 
s[[ c(3,2,3,1) ]]

s[[c("red", "green")]]
s$red

# expression based (partial) matching of names with single brackets
s["re"]
s["^re"]

# not with double brackets
# s[["re"]]


### SpatVector

v <- vect(system.file("ex/lux.shp", package="terra"))
v[2:3,]
v[1:2, 2:3]

subset(v, NAME_1 == "Diekirch", c(NAME_1, NAME_2))
}

\keyword{ spatial }

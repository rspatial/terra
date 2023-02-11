
\name{replace_values}

\docType{methods}

\alias{[<-}
\alias{[<-,SpatRaster,ANY,ANY,ANY-method}
\alias{[<-,SpatRaster,ANY,ANY-method}
\alias{[<-,SpatExtent,numeric,missing-method}

\alias{[<-,SpatVector,ANY,ANY-method}
\alias{[<-,SpatVector,ANY,missing-method}
\alias{[<-,SpatVector,missing,ANY-method}

\title{Replace values of a SpatRaster}

\description{
Replace values of a SpatRaster. These are convenience functions for smaller objects only. For larger rasters see \code{link{classify}} or \code{\link{subst}}
} 

\usage{
\S4method{[}{SpatRaster,ANY,ANY,ANY}(x, i, j, k) <- value

\S4method{[}{SpatVector,ANY,ANY}(x, i, j) <- value

\S4method{[}{SpatExtent,numeric,missing}(x, i, j) <- value
}

\arguments{
  \item{x}{SpatRaster}
  \item{i}{row numbers. numeric, logical, or missing for all rows. Can also be a SpatRaster or SpatVector}
  \item{j}{column numbers. numeric, logical or missing for all columns}
  \item{k}{layer number. numeric, logical or missing for all layers}
  \item{value}{numeric, matrix, or data.frame}
}
 
\value{SpatRaster}

\seealso{\code{\link{classify}, \link{subst}, \link{set.values}, \link{values}, \link{[[<-}}}

\examples{
## SpatRaster
r <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
r[] <- 1:25
r[1,] <- 5
r[,2] <- 10
r[r>10] <- NA

## SpatVector
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
v[2,2] <- "hello"
v[1,] <- v[10,]
v[,3] <- v[,1]
v[2, "NAME_2"] <- "terra"
head(v, 3)
}

\keyword{methods}
\keyword{spatial}


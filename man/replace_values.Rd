
\name{replace_values}

\docType{methods}

\alias{[<-,SpatRaster,ANY,ANY,ANY-method}

\title{Replace values of a SpatRaster}

\description{
Replace values of a SpatRaster. These are convenience functions for smaller objects only. For larger rasters see \code{link{classify}} or \code{\link{subst}}
} 

\usage{
\S4method{[}{SpatRaster,ANY,ANY,ANY}(x, i, j, k) <- value
}

\arguments{
  \item{x}{SpatRaster}
  \item{i}{row numbers. numeric, logical, or missing for all rows. Can also be a SpatRaster or SpatVector}
  \item{j}{column numbers. numeric, logical or missing for all columns}
  \item{k}{layer number. numeric, logical or missing for all layers}
  \item{value}{numeric, matrix, or data.frame}
}
 
\value{SpatRaster}

\seealso{\code{link{classify}, \link{subst}, \link{set.values}, \link{values}, \link{replace}}}

\examples{
r <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
r[] <- 1:25
r[1,] <- 5
r[,2] <- 10
r[r>10] <- NA
}

\keyword{methods}
\keyword{spatial}


\name{replace}

\docType{methods}

\alias{$<-,SpatRaster-method}
\alias{[,SpatRaster,logical,missing-method}
\alias{[<-,SpatRaster,logical,missing-method}
\alias{[<-,SpatRaster,missing,missing-method}
\alias{[<-,SpatRaster,missing,numeric-method}
\alias{[<-,SpatRaster,numeric,missing-method}
\alias{[<-,SpatRaster,numeric,numeric-method}
\alias{[<-,SpatRaster,SpatRaster,ANY-method}
\alias{[<-,SpatRaster,SpatVector,missing-method}
\alias{[[<-,SpatRaster,character,missing-method}
\alias{[[<-,SpatRaster,numeric,missing-method}

\title{Replace values of a SpatRaster}

\description{
Replace values of a SpatRaster. These are convenience functions for smaller objects only. For larger rasters see \code{link{classify}}
} 


\value{SpatRaster}

\seealso{\code{link{classify}, \link{values}, \link{replace}}}

\examples{
r <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
r[] <- 1:25
r[1,] <- 5
r[,2] <- 10
r[r>10] <- NA
}

\keyword{methods}
\keyword{spatial}


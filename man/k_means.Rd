\name{k_means}

\alias{k_means}

\alias{k_means,ANY-method}
\alias{k_means,SpatRaster-method}

\title{k_means}

\description{
Compute kmeans clusters for a SpatRaster
}

\usage{
\S4method{k_means}{SpatRaster}(x, centers, ..., filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{centers}{either the number of clusters, or a set of initial (distinct) cluster centres. If a number, a random set of (distinct) cells in \code{x} is chosen as the initial centres}   
  \item{...}{additional arguments passed to \code{\link[stats]{kmeans}}} 
  \item{maxcell}{positive integer. The size of the regular sample used if it is smaller than \code{ncell(x)}}
  \item{filename}{character. Output filename (ignored if \code{as.raster=FALSE}}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{ \code{\link[stats]{k_means}}
}

\examples{
f <- system.file("ex/logo.tif", package = "terra")
r <- rast(f)
km <- k_means(r, centers=5)
km
}

\keyword{spatial}


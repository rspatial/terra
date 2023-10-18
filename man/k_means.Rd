\name{k_means}

\alias{k_means}

\alias{k_means,ANY-method}
\alias{k_means,SpatRaster-method}

\title{k_means}

\description{
Compute kmeans clusters for a SpatRaster
}

\usage{
\S4method{k_means}{SpatRaster}(x, ..., filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{...}{additional arguments passed to \code{\link[stats]{k_means}}, such as the required argument \code{centers}} 
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


\name{princomp}

\alias{princomp}
\alias{princomp,SpatRaster-method}

\title{SpatRaster PCA}

\description{
Compute principal components for SpatRaster layers. 
The benefit of this method is that it uses all value to compute the principal components. 
Alternatively, you can use \code{stats::\link[stats]{princomp}} or \code{stats::\link[stats]{prcomp}} directly, perhaps using a sample. 
}

\usage{
\S4method{princomp}{SpatRaster}(x, cor=FALSE, nobs=FALSE)
}

\arguments{
  \item{x}{SpatRaster}
  \item{cor}{logical. If \code{FALSE}, the covariance matrix is used. Otherwise the correlation matrix is used (should only be used if there are no constant variables}   
  \item{nobs}{logical. Compute the number of observations? Ignored (and always \code{TRUE}) when \code{cor=TRUE}}
}

\value{
princomp object
}

\author{Based on a similar method by Benjamin Leutner}


\seealso{ \code{\link[stats]{princomp}}}

\examples{
f <- system.file("ex/logo.tif", package = "terra")
r <- rast(f)
pca <- princomp(r)
x <- predict(r, pca)

# use "index" to get a subset of the components
y <- predict(r, pca, index=1:2)

### use prcomp or princomp directly
pca <- princomp(r)

# may need to use sampling with a large raster 
sr <- spatSample(r, 100000, "regular")
pca <- prcomp(sr)
}

\keyword{spatial}


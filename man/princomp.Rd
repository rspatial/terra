\name{princomp}

\alias{princomp}
\alias{princomp,SpatRaster-method}

\title{SpatRaster PCA}

\description{
Compute principal components for SpatRaster layers. The benefit of this method is that it uses all values to compute the principal components, even for very large rasters. Alternatively, you can use \code{stats::\link[stats]{princomp}} or \code{stats::\link[stats]{prcomp}} directly, perhaps using a sample. 
}

\usage{
\S4method{princomp}{SpatRaster}(x, cor=FALSE, fix_sign=TRUE)
}

\arguments{
  \item{x}{SpatRaster}
  \item{cor}{logical. If \code{FALSE}, the covariance matrix is used. Otherwise the correlation matrix is used}   
  \item{fix_sign}{logical.  If \code{TRUE}, the signs of the loadings and scores are chosen so that the first element of each loading is non-negative}
}

\value{
princomp object
}

\author{Based on a similar method by Benjamin Leutner}

\note{
Use \code{mask(x, anyNA(x), maskvalue=TRUE)} if you want to assure that values from the same cells are used for all layers. Otherwise, the covariance matrix is computed for all cells in pairs of layers that are not \code{NA} in those two layers.
}

\seealso{ \code{\link[stats]{princomp}}}

\examples{
f <- system.file("ex/logo.tif", package = "terra")
r <- rast(f)
pca <- princomp(r)
x <- predict(r, pca)

# use "index" to get a subset of the components
p <- predict(r, pca, index=1:2)

### use princomp directly (or use prcomp instead)
pca <- princomp(r)
# may need to use sampling with a large raster 
sr <- spatSample(r, 100000, "regular")
pca <- prcomp(sr)
}

\keyword{spatial}


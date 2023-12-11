\name{princomp}

\alias{princomp}
\alias{princomp,SpatRaster-method}

\title{SpatRaster PCA}

\description{
Compute principal components for SpatRaster layers. The benefit of this method is that it uses all values to compute the principal components, even for very large rasters. This is because it computes the covariance matrix by processing the data in chunks, if necessary, using \code{\link{layerCor}}. The population covariance is used (not the sample, with \code{n-1} denominator, covariance). 

Alternatively, you can use \code{stats::\link[stats]{princomp}} or \code{stats::\link[stats]{prcomp}} with a data.frame of the raster values, using a sample for very large rasters (see the examples below). 
}

\usage{
\S4method{princomp}{SpatRaster}(x, cor=FALSE, fix_sign=TRUE, use="pairwise.complete.obs", maxcell=Inf)
}

\arguments{
  \item{x}{SpatRaster}
  \item{cor}{logical. If \code{FALSE}, the covariance matrix is used. Otherwise the correlation matrix is used}   
  \item{fix_sign}{logical.  If \code{TRUE}, the signs of the loadings and scores are chosen so that the first element of each loading is non-negative}  
  \item{use}{character. To decide how to handle missing values. This must be (an abbreviation of) one of the strings "everything", "complete.obs", "pairwise.complete.obs", or "masked.complete". With "pairwise.complete.obs", the covariance between a pair of layers is computed for all cells that are not \code{NA} in that pair. Therefore, it may be that the (number of) cells used varies between pairs. The benefit of this approach is that all available data is used. Use "complete.obs", if you want to only use the values from cells that are not \code{NA} in any of the layers. By using "masked.complete" you indicate that all layers have NA values in the same cells}
  \item{maxcell}{positive integer. The maximum number of cells to be used. If this is smaller than ncell(x), a regular sample of \code{x} is used}  
}

\value{
princomp object
}

\author{Alex Ilich and Robert Hijmans, based on a similar method by Benjamin Leutner}


\seealso{ \code{\link[stats]{princomp}}}

\examples{
f <- system.file("ex/logo.tif", package = "terra")
r <- rast(f)
pca <- princomp(r)
x <- predict(r, pca)

# use "index" to get a subset of the components
p <- predict(r, pca, index=1:2)

### use princomp directly
pca2 <- princomp(values(r),  fix_sign = TRUE)
p2 <- predict(r, pca2)

### may need to use sampling with a large raster 
### here with prcomp instead of princomp
sr <- spatSample(r, 100000, "regular")
pca3 <- prcomp(sr)
p3 <- predict(r, pca3)
}

\keyword{spatial}


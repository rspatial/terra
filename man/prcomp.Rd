\name{prcomp}

\alias{prcomp}
\alias{prcomp,SpatRaster-method}

\title{SpatRaster PCA with prcomp}

\description{
Compute principal components for SpatRaster layers. This method may be preferred to \code{\link{princomp}} for its greater numerical accuracy. However, it is slower and for very large rasters it can only be done with a sample. This may be good enough but see \code{\link{princomp}} if you want to use all values. Unlike \code{\link{princomp}}, in this method the sample variances are used with \code{n-1} as the denominator.  
}

\usage{
\S4method{prcomp}{SpatRaster}(x, retx=TRUE, center=TRUE, scale.=FALSE, 
		tol=NULL, rank.=NULL, maxcell=Inf)
}

\arguments{
  \item{x}{SpatRaster}
  \item{retx}{a logical value indicating whether the rotated variables should be returned}
  \item{center}{a logical value indicating whether the variables should be shifted to be zero centered. Alternately, a vector of length equal the number of columns of x can be supplied. The value is passed to \code{\link{scale}}}
  \item{scale.}{a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place. The default is FALSE for consistency with S, but in general scaling is advisable. Alternatively, a vector of length equal the number of columns of x can be supplied. The value is passed to \code{\link{scale}}}
  \item{tol}{a value indicating the magnitude below which components should be omitted. (Components are omitted if their standard deviations are less than or equal to tol times the standard deviation of the first component.) With the default null setting, no components are omitted (unless \code{rank.} is specified less than \code{min(dim(x))}). Other settings for \code{tol} could be \code{tol = 0} or \code{tol = sqrt(.Machine$double.eps)}, which would omit essentially constant components}
  \item{rank.}{optionally, a number specifying the maximal rank, i.e., maximal number of principal components to be used. Can be set as alternative or in addition to tol, useful notably when the desired rank is considerably smaller than the dimensions of the matrix}
  \item{maxcell}{positive integer. The maximum number of cells to be used. If this is smaller than ncell(x), a regular sample of \code{x} is used}  
}

\value{
prcomp object
}


\seealso{ \code{\link{princomp}}, \code{\link[stats]{prcomp}}}

\examples{
f <- system.file("ex/logo.tif", package = "terra")
r <- rast(f)
pca <- prcomp(r)
x <- predict(r, pca)

# use "index" to get a subset of the components
p <- predict(r, pca, index=1:2)

}

\keyword{spatial}


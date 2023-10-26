\name{layerCor}

\alias{layerCor}
\alias{layerCor,SpatRaster-method}


\title{Correlation and (weighted) covariance}


\description{
Compute correlation, (weighted) covariance, or similar summary statistics that compare the values of all pairs of the layers of a SpatRaster. 
}


\usage{
\S4method{layerCor}{SpatRaster}(x, fun, w, asSample=TRUE, na.rm=FALSE, maxcell=Inf, ...)
}


\arguments{
  \item{x}{SpatRaster}

  \item{fun}{character. The statistic to compute: either "cov" (covariance), "weighted.cov" (weighted covariance), or "pearson" (correlation coefficient) or your own function that takes two vectors as argument to compute a single number}

  \item{w}{SpatRaster with the weights  to compute the weighted covariance. It should have a single layer and the same geometry as \code{x}}

  \item{asSample}{logical. If \code{TRUE}, the statistic for a sample (denominator is \code{n-1}) is computed, rather than for the population (denominator is \code{n}). Only for the standard functions}
    
  \item{na.rm}{logical. Should missing values be removed?}
  
  \item{maxcell}{positive integer. The number of cells to be regularly sampled. Only used when \code{fun} is a function}

  \item{...}{additional arguments for \code{fun} (if it is a proper function)}
}


\value{
If \code{fun} is one of the three standard statistics, you get a list with two items: the correlation or (weighted) covariance matrix, and the (weighted) means.

If \code{fun} is a function, you get a matrix.
}

\references{
For the weighted covariance:
\itemize{
\item {Canty, M.J. and A.A. Nielsen, 2008. Automatic radiometric normalization of multitemporal satellite imagery with the iteratively re-weighted MAD transformation. Remote Sensing of Environment 112:1025-1036.}
\item {Nielsen, A.A., 2007. The regularized iteratively reweighted MAD method for change detection in multi- and hyperspectral data. IEEE Transactions on Image Processing 16(2):463-478.} 
}
}

\seealso{
 \code{\link{global}}, \code{\link{cov.wt}}, \code{\link{weighted.mean}}
}

\examples{
b <- rast(system.file("ex/logo.tif", package="terra"))   
layerCor(b, "pearson")

layerCor(b, "cov")

# weigh by column number
w <- init(b, fun="col")
layerCor(b, "weighted.cov", w=w)
}

	
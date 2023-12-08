\name{layerCor}

\alias{layerCor}
\alias{layerCor,SpatRaster-method}


\title{Correlation and (weighted) covariance}


\description{
Compute correlation, (weighted) covariance, or similar summary statistics that compare the values of all pairs of the layers of a SpatRaster. 
}


\usage{
\S4method{layerCor}{SpatRaster}(x, fun, w, asSample=TRUE, use="everything", maxcell=Inf, ...)
}


\arguments{
  \item{x}{SpatRaster}

  \item{fun}{character. The statistic to compute: either "cov" (covariance), "weighted.cov" (weighted covariance), or "cor" (pearson correlation coefficient) or your own function that takes two vectors as argument to compute a single number}

  \item{w}{SpatRaster with the weights  to compute the weighted covariance. It should have a single layer and the same geometry as \code{x}}

  \item{asSample}{logical. If \code{TRUE}, the statistic for a sample (denominator is \code{n-1}) is computed, rather than for the population (denominator is \code{n}). Only for the standard functions}
   
  \item{use}{character. To decide how to handle missing values. This must be (an abbreviation of) one of "everything", "complete.obs", "pairwise.complete.obs", "masked.complete". With "pairwise.complete.obs", the value for a pair of layers is computed for all cells that are not \code{NA} in that pair. Therefore, it may be that the (number of) cells used varies between pairs. The benefit of this approach is that all available data is used. Use "complete.obs", if you want to only use the values from cells that are not \code{NA} in any of the layers. By using "masked.complete" you indicate that all layers have NA values in the same cells}
    
  \item{maxcell}{positive integer. The maximum number of cells to be used. If this is smaller than ncell(x), a regular sample of \code{x} is used}

  \item{...}{additional arguments for \code{fun} (if it is a proper function)}
}


\value{
If \code{fun} is one of the three standard statistics, you get a list with two items: the correlation or (weighted) covariance matrix, and the (weighted) means. The means are also a matrix because they may depend on the combination of layers if different cells have NAs and are excluded from the computation. The rows of the mean matrix represent the layer whose (weighted) mean is being calculated and the columns represent the layer it is being paired with. Only cells with complete pair observations for both layers are used in the calculation of the (weighted) mean. The diagonal of the mean matrices are assigned NA.

If \code{fun} is a function, you get a single matrix.
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

	

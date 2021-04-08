\name{weighted.mean}

\alias{weighted.mean}
\alias{weighted.mean,SpatRaster,numeric-method}
\alias{weighted.mean,SpatRaster,SpatRaster-method}


\title{Weighted mean of layers}


\description{
Compute the weighted mean for each cell of the layers of a SpatRaster. The weights can be spatially variable or not.
}


\usage{
\S4method{weighted.mean}{SpatRaster,numeric}(x, w, na.rm=FALSE, filename="", ...)

\S4method{weighted.mean}{SpatRaster,SpatRaster}(x, w, na.rm=FALSE, filename="", ...)
}


\arguments{
  \item{x}{SpatRaster}

  \item{w}{A vector of weights (one number for each layer), or for spatially variable weights, a SpatRaster with weights (should have the same extent, resolution and number of layers as x)}

  \item{na.rm}{Logical. Should missing values be removed?}

  \item{filename}{character. Output filename}
  \item{...}{options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}


\seealso{
\code{\link{Summary-methods}}, \code{\link[stats]{weighted.mean}}
}


\examples{
b <- rast(system.file("ex/logo.tif", package="terra"))   

# give least weight to first layer, most to last layer
wm1 <- weighted.mean(b, w=1:3)

# spatially varying weights
# weigh by column number
w1 <- init(b, "col")

# weigh by row number
w2 <- init(b, "row")
w <- c(w1, w2, w2)

wm2 <- weighted.mean(b, w=w)
}


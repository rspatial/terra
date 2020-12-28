\name{pairs}

\docType{methods}

\alias{pairs}
\alias{pairs,SpatRaster-method}

\title{
Pairs plot (matrix of scatterplots)
}

\description{
Pair plots of layers in a SpatRaster. This is a wrapper around graphics function \code{\link[graphics]{pairs}}.
}

\usage{
\S4method{pairs}{SpatRaster}(x, hist=TRUE, cor=TRUE, use="pairwise.complete.obs", maxcells=100000, ...)
}

\arguments{
\item{x}{SpatRaster}
\item{hist}{logical. If TRUE a histogram of the values is shown on the diagonal}
\item{cor}{logical. If TRUE the correlation coefficient is shown in the upper panels}
\item{use}{argument passed to the \code{\link[stats]{cor}} function}
\item{maxcells}{integer. Number of pixels to sample from each layer of a large SpatRaster}
\item{...}{additional arguments (graphical parameters)}
}

\seealso{
\code{\link{boxplot}, \link{hist}}
}

\examples{
r <-rast(system.file("ex/elev.tif", package="terra"))
s <- c(r, 1/r, sqrt(r))
names(s) <- c("elevation", "inverse", "sqrt") 
pairs(s)

# to make indvidual histograms:
hist(r)
# or scatter plots:
plot(s[[1]], s[[2]])
}


\keyword{spatial}


\name{terra-package}
\alias{terra-package}

\alias{terra}

\docType{package}

\title{terra}

\description{
This package is conceived as a replacement of the raster package. It has a very similar, but simpler, interface, and it is faster.

An important difference is that there is a single main class, "SpatRaster", as opposed to the three classes in the raster package (RasterLayer, RasterStack, RasterBrick). 

To speed up computations, all native computations are done in C++, exposed via a Rcpp module. Methods that heavily interact with native R functions (e.g. "predict") are written in R, calling C++ functions for reading/writing.
}


\keyword{ package }
\keyword{ spatial }


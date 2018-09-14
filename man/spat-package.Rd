\name{terra-package}
\alias{terra-package}

\alias{terra}

\docType{package}

\title{terra}

\description{
This package is conceived as a replacement of the raster package. It has a very similar, but simpler, interface, and it is faster.
 
There is a single main class, "SpatRaster", as opposed to the three classes in the raster package (RasterLayer, Stack, Brick). 

To speed up computations, all native computations are done in C++. 

Methods that heavily interact with native R functions (e.g. "predict") will be written in R, calling C++ functions for reading/writing.

Classes, methods and properties are exposed via a Rcpp module.

There is one main S4 class "SpatRaster" that holds a C++ reference class to a C++ SpatRaster object.

The normal end-user won't directly use the reference class.

A new design rule each algorithm is implemented only once... In the raster package many algos are implemented twice, once for in memory, once for chunk wise to disk computations for large files. 

Note that the C++ code uses 0 based indexing for row, col, layer, and cell numbers. 

Raster data are three-dimensional (and in the future there could be support for more). But internally they are stored as one dimensional vectors. But I may use armadillo matrix and cube classes within methods.

}


\keyword{ package }
\keyword{ spatial }


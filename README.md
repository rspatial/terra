# terra

This package is conceived as a replacement of the raster package. It has a very similar, but simpler, interface, and it is faster.

To speed up computations, all native computations are done in C++. Classes, methods and properties are exposed via a Rcpp module. The normal end-user won't directly use the reference classes. There is a single main class for raster data, "SpatRaster". This S4 class holds a C++ reference class to a C++ SpatRaster object.

It is still in a early stage of development. 

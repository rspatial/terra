# terra

This R package is a replacement of the raster package. It has a very similar, but has a simpler interface, and it is faster.

All native computations are done in C++. Classes, methods and properties are exposed via a Rcpp module. The R side has three S4 classes (SpatRaster, SpatVector and SpatExtent) that hold a reference to a C++ object. 


# terra

This R package is a replacement of the raster package. It has a very similar, but has a simpler interface, and it is faster.

All native computations are done in C++. Classes, methods and properties are exposed via a Rcpp module. The R side has three S4 classes (SpatRaster, SpatVector and SpatExtent) that hold a reference to a C++ object. 

## Installation


### Windows

If you are on Windows, you need to first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) to get a C++ compiler that R can use. After that, you need the first install development version of "raster" for "terra" to work. 

Then, in R, install the packages.

```
library(devtools)
devtools::install_github("rspatial/raster")
devtools::install_github("rspatial/terra")
```

### Mac - OSX

The libraries GDAL (>= 2.0.0), GEOS (>= 3.3.0) and Proj.4 (>= 4.8.0) are required (as for [https://github.com/r-spatial/sf](sf) )

With Homebrew you can do

```
brew unlink gdal
brew tap osgeo/osgeo4mac && brew tap --repair
brew install proj
brew install geos
brew install gdal2 --with-armadillo --with-complete --with-libkml --with-unsupported
brew link --force gdal2
```

And now, in R, install the packages.
```
library(devtools)
devtools::install_github("rspatial/raster")
devtools::install_github("rspatial/terra")
```

### Linux

The libraries GDAL (>= 2.0.0), GEOS (>= 3.3.0) and Proj.4 (>= 4.8.0) are required (as for [https://github.com/r-spatial/sf](sf) )


To install these on Ubuntu you can do:
```
sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
sudo apt-get update
sudo apt-get install libgdal-dev libgeos-dev libproj-dev 
```

And now, in R, install the packages.
```
library(devtools)
devtools::install_github("rspatial/raster")
devtools::install_github("rspatial/terra")
```

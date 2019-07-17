# terra

This R package is a replacement of the [raster](https://github.com/rspatial/raster) package.
It has a very similar, but simpler, interface, and it is much faster.

All native computations are done in C++. 
Classes, methods and properties are exposed via a Rcpp module.
The R side has two main S4 classes (SpatRaster and SpatVector) that represent spatial data. These classes have only slot, a reference to a C++ object. They are used to provide a "normal" "S4" R user-interface as in the raster package.

The first (alpha) release is expected by July 2019.

## Installation

You need to have the current version of `raster` from CRAN (>= 2.9-22) or install the [development version of "raster"](https://github.com/rspatial/raster).

### Windows

If you are on Windows, you need to first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) to get a C++ compiler that R can use.
After that, you need the first install development version of "raster" for "terra" to work. 

Then, in R, install the packages.

```
library(devtools)
install.packages("raster")
#devtools::install_github("rspatial/raster")
devtools::install_github("rspatial/terra")
```

### Mac - OSX

The libraries GDAL (>= 2.0.0), GEOS (>= 3.3.0) and Proj.4 (>= 4.8.0) are required (as for [sf](https://github.com/r-spatial/sf))

With Homebrew you can do:

```
brew install gdal
```

And now, in R, install the packages.
```
library(devtools)
install.packages("raster")
#devtools::install_github("rspatial/raster")
devtools::install_github("rspatial/terra")
```

### Linux

The libraries GDAL (>= 2.0.0), GEOS (>= 3.3.0) and Proj.4 (>= 4.8.0) are required (as for [sf](https://github.com/r-spatial/sf))


To install these on Ubuntu you can do:
```
sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
sudo apt-get update
sudo apt-get install libgdal-dev libgeos-dev libgeos++-dev libproj-dev 
```

And now, in R, install the packages.
```
library(devtools)
install.packages("raster")
#devtools::install_github("rspatial/raster")
devtools::install_github("rspatial/terra")
```

See the sf [instructions](https://github.com/r-spatial/sf) for installation on other linux systems.


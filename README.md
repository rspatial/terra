# terra

[![Travis build
status](https://travis-ci.org/rspatial/terra.svg?branch=master)](https://travis-ci.org/rspatial/terra)

`terra` is a R package is a replacement of the [raster](https://github.com/rspatial/raster) package.
It has a very similar, but simpler, interface, and it is much faster.

All native computations are done in C++.  Classes, methods and properties are exposed via a Rcpp module.
The R side has two main S4 classes (SpatRaster and SpatVector) that represent spatial data. These classes have only slot, a reference to a C++ object. They are used to provide a "normal" "S4" R user-interface as in the raster package.


## Installation

`terra` is available from CRAN, so you can use `install.packages("terra")`. 

See below for instructions on installing the *development version*

### Windows

If you are on Windows, you need to first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) to get a C++ compiler that R can use. After that, you need the first install development version of "raster" for "terra" to work. 

Then, in R, install the packages.

```
remotes::install_github("rspatial/terra")
```

### Mac - OSX

The libraries GDAL (>= 2.0.0), GEOS (>= 3.3.0) and Proj.4 (>= 4.8.0) are required (as for [sf](https://github.com/r-spatial/sf))

Install Homebrew, and then, with Homebrew you can do, from the terminal:

```
brew install pkg-config
brew install gdal
```

And now, in R, install the package
```
remotes::install_github("rspatial/terra")
```

If you get `configure: error: libproj not found in standard or given locations.` Then you need to clone the repo and install like this

```
remotes::install_github("rspatial/terra", configure-args="--with-proj-lib=/usr/local/lib/")
```

If you get error `unknown type name 'uuid_t'` you first need to install Rcpp > 1.0.4 (currently, the development version from github)

```
remotes::install_github("RcppCore/Rcpp.git") 
```


### Linux

The libraries GDAL (>= 2.0.0), GEOS (>= 3.3.0) and Proj.4 (>= 4.8.0) are required (as for [sf](https://github.com/r-spatial/sf))


To install these on Ubuntu you can do:
```
sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
sudo apt-get update
sudo apt-get install libgdal-dev libgeos-dev libproj-dev 
```

And now, in R, install the packages.
```
remotes::install_github("rspatial/terra")
```

See the sf [instructions](https://github.com/r-spatial/sf) for installation on other linux systems --- and for possible updates/improvements on the above instructions.


# terra

[![Travis build
status](https://travis-ci.org/rspatial/terra.svg?branch=master)](https://travis-ci.org/rspatial/terra)
[![CRAN
status](https://www.r-pkg.org/badges/version/terra)](https://cran.r-project.org/package=terra)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/terra)](http://www.r-pkg.org/pkg/terra)


`terra` is an R package that replaces [raster](https://github.com/rspatial/raster).
It has a very similar interface, but it is simpler and much faster. The first (beta-) release was on 20 March 2020.

`terra` is written in C++.  Classes, methods and properties are exposed via a Rcpp module. The R side has two main classes (SpatRaster and SpatVector) that represent spatial data. These classes are used to provide a standard R user-interface. There are tutorials at [rspatial.org/terra](https://rspatial.org/terra/index.html).


## Installation

`terra` is available from CRAN, so you can use `install.packages("terra")`. 

See below for instructions on installing the *development version*

### All OS

First intstall the packages that terra depends on 

```
install.packages(c("raster", "Rcpp"))
```

### Windows

If you are on Windows, you need to first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) to get a C++ compiler that R can use. 

Then, in R, install the package.

```
remotes::install_github("rspatial/terra")
```


### MacOS

First install gdal and proj with homebrew

```
brew install pkg-config
brew install gdal
```

Installation may require additional configuration, that you can pass on with `devtools::install_github`

```
library(devtools)
devtools::install_github("rspatial/terra", configure.args = "--with-proj-lib=/usr/local/lib/")
```

### Linux

The libraries GDAL (>= 3.0.0), GEOS (>= 3.3.0) and Proj.4 (>= 6.0.0) are required 


To install these on Ubuntu version 18.04 (Bionic) you can do:
```
sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
sudo apt-get update
sudo apt-get install libgdal-dev libgeos-dev libproj-dev 
```

And now, in R, install the package
```
remotes::install_github("rspatial/terra")
```

See the `sf` [instructions](https://github.com/r-spatial/sf) for installation on other linux systems --- and for possible updates/improvements on the above instructions. But note that `terra` depends on on more recent versions of GDAL and PROJ libraries 


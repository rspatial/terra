# terra

[![rcmdcheck](https://github.com/rspatial/terra/actions/workflows/rcmdcheck.yml/badge.svg)](https://github.com/rspatial/terra/actions/workflows/rcmdcheck.yml)
[![CRAN
status](https://www.r-pkg.org/badges/version/terra)](https://cran.r-project.org/package=terra)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/terra)](http://www.r-pkg.org/pkg/terra)

<p align="right"> 
  <img align="right" width="350" height="350" src="https://github.com/rspatial/terra/raw/master/logo.png">
logo by Zane Dax</p>

`terra` is an R package for spatial analysis. There are tutorials at [rspatial.org/terra](https://rspatial.org/terra/index.html). 

[stackoverflow](https://stackoverflow.com/) is the best place to ask questions if you get stuck. Make sure to include a [simple reproducible example](https://stackoverflow.com/questions/5963269/how-to-make-a-great-r-reproducible-example). But if you think you have found a bug, please file an [issue](https://github.com/rspatial/terra/issues).

`terra` replaces the [raster](https://github.com/rspatial/raster) package. The interfaces of `terra` and `raster` are similar, but `terra` is simpler, faster and can do more. 


## Installation

`terra` is available from CRAN, so you can use `install.packages("terra")` to get the current *released version*

The easiest way to use the *development version* on Windows or MacOS, is to install it from the [R-universe](https://r-universe.dev/organizations/), like this:

```
install.packages('terra', repos='https://rspatial.r-universe.dev')
```

Please note that, on MacOS, the development version does not support the HDF4 file format, which is used for some NASA satellite data products.

### From source-code

To install from source-code, first install the packages that terra depends on: 

```
install.packages(c("raster", "Rcpp"))
```

And then continue based on the OS you are using. 

#### Windows

On Windows, you need to first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) to get a C++ compiler that R can use. 

Then, in R, install the package.

```
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
remotes::install_github("rspatial/terra")
```

#### MacOS

On OSX, first install gdal and proj with homebrew

```
brew install pkg-config
brew install gdal
```
Followed by

```
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
remotes::install_github("rspatial/terra")
```

This should work on **Catalina** and **Big Sur**


#### Linux

The GDAL (>= 2.2.0), GEOS (>= 3.3.0) and PROJ (>= 6.0.0) libraries are required 


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


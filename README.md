# terra

[![Travis build
status](https://travis-ci.org/rspatial/terra.svg?branch=master)](https://travis-ci.org/rspatial/terra)
[![CRAN
status](https://www.r-pkg.org/badges/version/terra)](https://cran.r-project.org/package=terra)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/terra)](http://www.r-pkg.org/pkg/terra)


`terra` is an R package that replaces [raster](https://github.com/rspatial/raster).
It has a very similar interface, but it is simpler, much faster and can do more. It was first released on 20 March 2020. There are tutorials at [rspatial.org/terra](https://rspatial.org/terra/index.html).


## Installation

### CRAN 

`terra` is available from CRAN, so you can use `install.packages("terra")` to get the current *released version*

### R-Universe

The easiest  way to install the *development version* is to get it from r-universe like this

```
options(repos = c(
    rspatial = 'https://rspatial.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))

install.packages('terra')
```

### Self compilation

For a more bare bones install, first install the packages that terra depends on 

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


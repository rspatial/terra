# terra

[![CRAN
status](https://www.r-pkg.org/badges/version/terra)](https://cran.r-project.org/package=terra)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/terra)](http://www.r-pkg.org/pkg/terra)

`terra` is an R package for spatial analysis. There are tutorials at [rspatial.org/terra](https://rspatial.org/terra/index.html). 

[stackoverflow](https://stackoverflow.com/) is the best place to ask questions if you get stuck. Make sure to include a [simple reproducible example](https://stackoverflow.com/questions/5963269/how-to-make-a-great-r-reproducible-example). But if you think you have found a bug, please file an [issue](https://github.com/rspatial/raster/issues).

`terra` replaces [raster](https://github.com/rspatial/raster). The interfaces of `terra` and `raster` are similar, but `terra` it is simpler, much faster and it can do more. 

## Installation

### CRAN 

`terra` is available from CRAN, so you can use `install.packages("terra")` to get the current *released version*

The easiest way to install the *development version* on Windows or MacOS, is to get it from the [R-universe](https://r-universe.dev/organizations/) repository, like this:

```
install.packages('terra', repos='https://rspatial.r-universe.dev')
```

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


# terra <img align="right" width="250" src="man/figures/logo.png">

<p align="right"; style="font-size:11px">logo by Zane Dax</p>

[![rcmdcheck](https://github.com/rspatial/terra/actions/workflows/rcmdcheck.yml/badge.svg)](https://github.com/rspatial/terra/actions/workflows/rcmdcheck.yml)
[![CRAN
status](https://www.r-pkg.org/badges/version/terra)](https://cran.r-project.org/package=terra)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/terra)](http://www.r-pkg.org/pkg/terra)

`terra` is an *R* package for spatial data analysis. There are tutorials at [rspatial.org/terra](https://rspatial.org/terra/index.html). 

[Stack Overflow](https://stackoverflow.com/questions/tagged/terra) is the best place to ask questions if you get stuck. Make sure to include a [simple reproducible example](https://stackoverflow.com/questions/5963269/how-to-make-a-great-r-reproducible-example). But if you think you have found a bug, please file an [issue](https://github.com/rspatial/terra/issues).

`terra` replaces the [raster](https://github.com/rspatial/raster) package. The interfaces of `terra` and `raster` are similar, but `terra` is simpler, faster and can do more. 


## Installation

`terra` is available from CRAN, so you can use `install.packages("terra")` to get the current *released version*.

The easiest way to use the *development version* on Windows or MacOS, is to install it from the [R-universe](https://r-universe.dev/organizations/), like this:


```
install.packages('terra', repos='https://rspatial.r-universe.dev')
```


### From source-code

To install from source-code, first install the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) package that terra depends on: 

```
install.packages("Rcpp")
```

And then continue based on the OS you are using. 

#### Windows

On Windows, you need to first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) to get a C++ compiler that R can use. You need a recent version of Rtools42 (rtools42-5355-5357).

Then, in R, install the package.

```
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
remotes::install_github("rspatial/terra")
```

#### macOS

On macOS, you can use [MacPorts](https://www.macports.org/) or [Homebrew](https://brew.sh/).

With MacPorts you can do 

```
sudo port install R-terra
```

With Homebrew, you need to first install GDAL:

```
brew install pkg-config
brew install gdal
```

Followed by (note the additional configuration argument needed for Homebrew)

```
remotes::install_github("rspatial/terra", configure.args = "--with-proj-lib=$(brew --prefix)/lib/")
```

To install the CRAN version from source you would do

```
install.packages("terra", type = "source", configure.args = "--with-proj-lib=$(brew --prefix)/lib/")
```

#### Linux

The *easy* way to install terra on Ubuntu is with [r2u](https://eddelbuettel.github.io/r2u/).

The harder way: 

Install the system requirements GDAL (>= 2.2.3), GEOS (>= 3.4.0), PROJ (>= 4.9.3), sqlite3.


```
sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
sudo apt-get update
sudo apt-get install libgdal-dev libgeos-dev libproj-dev 
```

And now, in R, install the package
```
remotes::install_github("rspatial/terra")
```

See the `sf` [instructions](https://github.com/r-spatial/sf) for installation on other linux systems --- and for possible updates/improvements on the above instructions.


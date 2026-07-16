# Path-Based Nondisperive Flow Direction

Computes non-disperive flow/draiage direction according to path-based
methods over grid-based digital elevation models.

## Usage

``` r
# S4 method for class 'SpatRaster'
flowdirD8ltd(x,lambda=0.5,deviation_type=c("ltd","lad"),max_iters=10^6,filename="", ...) 
# S4 method for class 'SpatRaster'
flowdirD8lad(x,lambda=0.5,deviation_type=c("lad","ltd"),max_iters=10^6,filename="",...)
```

## Arguments

- x:

  SpatRaster with elevation/terrain model , see
  [`terrain`](https://rspatial.github.io/terra/reference/terrain.md).

- lambda:

  Parameter between 0 and 1.

- deviation_type:

  Character string. Default is the first element of `=c("ltd","lad")`.
  If `"ltd"` (default) flow direction dispersion is deteceted with LTD
  criterion, if `"lad"` flow direction dispersion is deteceted with LTD
  criterion. See Orlandini et al,2003 for further details.

- max_iters:

  maximum iterations for drainage path starting points detection

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## Details

The algorithm is an adaptation of the one proposed by Li et al, 2021 and
Orlandini et al,2003 This function is experimental and under
development: results are to be verified.

## See also

[`terrain`](https://rspatial.github.io/terra/reference/terrain.md),[`watershed`](https://rspatial.github.io/terra/reference/watershed.md),[`flowAccumulation`](https://rspatial.github.io/terra/reference/flowAccumulation.md)

## Author

Zhenya Li at al, Stefano Orlandini at al. , Emanuele Cordano (R/C/C++
code implemetation)

## References

Orlandini, S., G. Moretti, M. Franchini, B. Aldighieri, and B. Testa
(2003), Path-based methods for the determination of nondispersive
drainage directions in grid-based digital elevation models, Water
Resour. Res., 39, 1144, doi:10.1029/2002WR001639, 6.
<https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2002WR001639>

Li, Z., Shi, P., Yang, T., Wang, C., Yong, B., & Song, Y. (2022). An
improved D8-LTD for the extraction of total contributing area (TCA) by
adopting the strategies of path independency and local dispersion. Water
Resources Research, 58, e2021WR030948.
https://doi.org/10.1029/2021WR030948
<https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021WR030948>

Useful presentation:
<http://www.idrologia.unimore.it/orlandini/web-archive/seminars/nyc-2008-2.pdf>

## Examples

``` r
elev1 <- array(NA,c(9,9))
elev2 <- elev1
dx <- 1
dy <- 1 
for (r in 1:nrow(elev1)) {
  y <- (r-5)*dx
  for (c in 1:ncol(elev1)) {
    
    x <- (c-5)*dy
    elev1[r,c] <- 5*(x^2+y^2)
    elev2[r,c] <- 10+5*(abs(x))-0.001*y 
  }
} 


## Elevation Raster Maps
elev1 <- rast(elev1)
elev2 <- rast(elev2)

t(array(elev1[],rev(dim(elev1)[1:2])))
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]  160  125  100   85   80   85  100  125  160
#>  [2,]  125   90   65   50   45   50   65   90  125
#>  [3,]  100   65   40   25   20   25   40   65  100
#>  [4,]   85   50   25   10    5   10   25   50   85
#>  [5,]   80   45   20    5    0    5   20   45   80
#>  [6,]   85   50   25   10    5   10   25   50   85
#>  [7,]  100   65   40   25   20   25   40   65  100
#>  [8,]  125   90   65   50   45   50   65   90  125
#>  [9,]  160  125  100   85   80   85  100  125  160
t(array(elev2[],rev(dim(elev2)[1:2])))
#>         [,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]   [,8]   [,9]
#>  [1,] 30.004 25.004 20.004 15.004 10.004 15.004 20.004 25.004 30.004
#>  [2,] 30.003 25.003 20.003 15.003 10.003 15.003 20.003 25.003 30.003
#>  [3,] 30.002 25.002 20.002 15.002 10.002 15.002 20.002 25.002 30.002
#>  [4,] 30.001 25.001 20.001 15.001 10.001 15.001 20.001 25.001 30.001
#>  [5,] 30.000 25.000 20.000 15.000 10.000 15.000 20.000 25.000 30.000
#>  [6,] 29.999 24.999 19.999 14.999  9.999 14.999 19.999 24.999 29.999
#>  [7,] 29.998 24.998 19.998 14.998  9.998 14.998 19.998 24.998 29.998
#>  [8,] 29.997 24.997 19.997 14.997  9.997 14.997 19.997 24.997 29.997
#>  [9,] 29.996 24.996 19.996 14.996  9.996 14.996 19.996 24.996 29.996

plot(elev1)

plot(elev2)


## Flow direction raster
flowdir1l1<- flowdirD8ltd(elev1,lambda=1)
flowdir2l1<- flowdirD8ltd(elev2,lambda=1)

t(array(flowdir1l1[],rev(dim(flowdir1l1)[1:2])))
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    2    2    4    4    4    4    4    8    8
#>  [2,]    2    2    2    4    4    4    8    8    8
#>  [3,]    1    2    2    4    4    4    8    8   16
#>  [4,]    1    1    1    2    4    8   16   16   16
#>  [5,]    1    1    1    1    0   16   16   16   16
#>  [6,]    1    1    1  128   64   32   16   16   16
#>  [7,]    1  128  128   64   64   64   32   32   16
#>  [8,]  128  128  128   64   64   64   32   32   32
#>  [9,]  128  128   64   64   64   64   64   32   32
t(array(flowdir2l1[],rev(dim(flowdir2l1)[1:2])))
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    1    1    1    1    4   16   16   16   16
#>  [2,]    1    1    1    1    4   16   16   16   16
#>  [3,]    1    1    1    1    4   16   16   16   16
#>  [4,]    1    1    1    1    4   16   16   16   16
#>  [5,]    1    1    1    1    4   16   16   16   16
#>  [6,]    1    1    1    1    4   16   16   16   16
#>  [7,]    1    1    1    1    4   16   16   16   16
#>  [8,]    1    1    1    1    4   16   16   16   16
#>  [9,]    1    1    1    1    0   16   16   16   16


flowdir1ldf<- flowdirD8ltd(elev1,lambda=0.5)
flowdir2ldf<- flowdirD8ltd(elev2,lambda=0.5)

t(array(flowdir1ldf[],rev(dim(flowdir1l1)[1:2])))
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    2    2    4    4    4    4    4    8    8
#>  [2,]    2    2    4    4    4    4    4    8    8
#>  [3,]    1    1    2    4    4    4    8   16   16
#>  [4,]    1    1    1    2    4    8   16   16   16
#>  [5,]    1    1    1    1    0   16   16   16   16
#>  [6,]    1    1    1  128   64   32   16   16   16
#>  [7,]    1    1  128   64   64   64   32   16   16
#>  [8,]  128  128   64   64   64   64   64   32   32
#>  [9,]  128  128   64   64   64   64   64   32   32
t(array(flowdir2ldf[],rev(dim(flowdir2l1)[1:2])))
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    1    1    1    1    4   16   16   16   16
#>  [2,]    1    1    1    1    4   16   16   16   16
#>  [3,]    1    1    1    1    4   16   16   16   16
#>  [4,]    1    1    1    1    4   16   16   16   16
#>  [5,]    1    1    1    1    4   16   16   16   16
#>  [6,]    1    1    1    1    4   16   16   16   16
#>  [7,]    1    1    1    1    4   16   16   16   16
#>  [8,]    1    1    1    1    4   16   16   16   16
#>  [9,]    1    1    1    1    0   16   16   16   16




flowdir1l0<- flowdirD8ltd(elev1,lambda=0)
flowdir2l0<- flowdirD8ltd(elev2,lambda=0)

t(array(flowdir1l0[],rev(dim(flowdir1l1)[1:2])))
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    2    2    4    4    4    4    4    8    8
#>  [2,]    2    2    2    4    4    4    8    8    8
#>  [3,]    1    2    2    4    4    4    8    8   16
#>  [4,]    1    1    1    2    4    8   16   16   16
#>  [5,]    1    1    1    1    0   16   16   16   16
#>  [6,]    1    1    1  128   64   32   16   16   16
#>  [7,]    1  128  128   64   64   64   32   32   16
#>  [8,]  128  128  128   64   64   64   32   32   32
#>  [9,]  128  128   64   64   64   64   64   32   32
t(array(flowdir2l0[],rev(dim(flowdir2l1)[1:2])))
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    1    1    1    1    4   16   16   16   16
#>  [2,]    1    1    1    1    4   16   16   16   16
#>  [3,]    1    1    1    1    4   16   16   16   16
#>  [4,]    1    1    1    1    4   16   16   16   16
#>  [5,]    1    1    1    1    4   16   16   16   16
#>  [6,]    1    1    1    1    4   16   16   16   16
#>  [7,]    1    1    1    1    4   16   16   16   16
#>  [8,]    1    1    1    1    4   16   16   16   16
#>  [9,]    1    1    1    1    0   16   16   16   16


## Flow direction raster
flowdir1<- terrain(elev1,v="flowdir")
flowdir2<- terrain(elev2,v="flowdir")


## Cone Geometry using calculus R package 


library(calculus)
#> 
#> Attaching package: ‘calculus’
#> The following object is masked from ‘package:terra’:
#> 
#>     wrap
library(stringr)

dx <- 1 #2.5
dy <- 1 #2.5
xmin <- -10
ymin <- -10

xmax <- 10
ymax <- 10


x <- seq(from=xmin,to=xmax,by=dx)
y <- seq(from=ymin,to=ymax,by=dy)


####
vars <- list()
vars$x <- rep(x,times=length(y))
vars$y <- rep(y,each=length(x))
vars <- as.data.frame(vars)
### 
elev_f <- "(x^2+y^2)^(1/2)" ## cone ##ok
###

vars$elev <- evaluate(elev_f,var=vars[,c("x","y")])
vars$elev[which(vars$x %in% range(vars$x))] <- NA 
vars$elev[which(vars$y %in% range(vars$y))] <- NA 

## Flow Angle 

vgrad <- matrix(gradient(elev_f,var=vars[,c("x","y")],accuracy=8),ncol=2)
vars$flow_angle <- atan2(y=-vgrad[,2],x=-vgrad[,1])+vars$elev*0


####


rr <- rast(vars)
rr$flow_dir <- terrain(rr$elev,"flowdir")
rr$flow_dirltdl0 <- flowdirD8ltd(rr$elev,lambda=0)
rr$flow_dirltdlm <- flowdirD8ltd(rr$elev,lambda=0.5)
rr$flow_dirltdl1 <- flowdirD8ltd(rr$elev,lambda=1)

plot(rr$elev)
arrows_on_rast(rr$flow_dir, unit="flowdir",col="black")

#> NULL
plot(rr$elev)
arrows_on_rast(rr$flow_dirltdl0, unit="flowdir",col="blue")

#> NULL
plot(rr$elev)
arrows_on_rast(rr$flow_dirltdlm, unit="flowdir",col="black")

#> NULL
plot(rr$elev)
arrows_on_rast(rr$flow_dirltdl1, unit="flowdir",col="black")

#> NULL

## V-Shape Watershed  Geometry using calculus R package 

library(calculus)
library(stringr)

dx <- 1 #2.5
dy <- 1 #2.5
xmin <- -10
ymin <- -10

xmax <- 10
ymax <- 10


x <- seq(from=xmin,to=xmax,by=dx)
y <- seq(from=ymin,to=ymax,by=dy)


####
vars <- list()
vars$x <- rep(x,times=length(y))
vars$y <- rep(y,each=length(x))
vars <- as.data.frame(vars)
### 
elev_f <- "x*tanh(100*x)+0.2*y" ## v-shape basin
###

vars$elev <- evaluate(elev_f,var=vars[,c("x","y")])
vars$elev[which(vars$x %in% range(vars$x))] <- NA 
vars$elev[which(vars$y %in% range(vars$y))] <- NA 

## Flow Angle 

vgrad <- matrix(gradient(elev_f,var=vars[,c("x","y")],accuracy=8),ncol=2)
vars$flow_angle <- atan2(y=-vgrad[,2],x=-vgrad[,1])+vars$elev*0


####


rr <- rast(vars)
rr$flow_dir <- terrain(rr$elev,"flowdir")
rr$flow_dirladl0 <- flowdirD8ltd(rr$elev,lambda=0,deviation_type="lad") 
##or flowdirD8lad(rr$elev,lambda=0)
rr$flow_dirladlm <- flowdirD8ltd(rr$elev,lambda=0.5,deviation_type="lad") 
##or flowdirD8lad(rr$elev,lambda=0.5)
rr$flow_dirladl1 <- flowdirD8ltd(rr$elev,lambda=1,deviation_type="lad") 
##or flowdirD8lad(rr$elev,lambda=1)
rr$flow_dirltdl0 <- flowdirD8ltd(rr$elev,lambda=0,deviation_type="ltd")
rr$flow_dirltdlm <- flowdirD8ltd(rr$elev,lambda=0.5,deviation_type="ltd")
rr$flow_dirltdl1 <- flowdirD8ltd(rr$elev,lambda=1,deviation_type="ltd")

plot(rr$elev)
arrows_on_rast(rr$flow_dir, unit="flowdir",col="black")

#> NULL
plot(rr$elev)
arrows_on_rast(rr$flow_dirladl0, unit="flowdir",col="blue")

#> NULL
plot(rr$elev)
arrows_on_rast(rr$flow_dirladlm, unit="flowdir",col="black")

#> NULL
plot(rr$elev)
arrows_on_rast(rr$flow_dirladl1, unit="flowdir",col="black")

#> NULL

plot(rr$elev)
arrows_on_rast(rr$flow_dirltdl0, unit="flowdir",col="blue")

#> NULL
plot(rr$elev)
arrows_on_rast(rr$flow_dirltdlm, unit="flowdir",col="black")

#> NULL
plot(rr$elev)
arrows_on_rast(rr$flow_dirltdl1, unit="flowdir",col="black")

#> NULL



## Planar Hillslope  Geometry using calculus R package 

library(calculus)
library(stringr)

dx <- 1 #2.5
dy <- 1 #2.5
xmin <- -10
ymin <- -10

xmax <- 10
ymax <- 10


x <- seq(from=xmin,to=xmax,by=dx)
y <- seq(from=ymin,to=ymax,by=dy)


####
vars <- list()
vars$x <- rep(x,times=length(y))
vars$y <- rep(y,each=length(x))
vars <- as.data.frame(vars)
### 
elev_f <- "x+0.5*y" ## planar hillslope basin , try also with elev_f <- "x+0.505*y"

###

vars$elev <- evaluate(elev_f,var=vars[,c("x","y")])
vars$elev[which(vars$x %in% range(vars$x))] <- NA 
vars$elev[which(vars$y %in% range(vars$y))] <- NA 

## Flow Angle 

vgrad <- matrix(gradient(elev_f,var=vars[,c("x","y")],accuracy=8),ncol=2)
vars$flow_angle <- atan2(y=-vgrad[,2],x=-vgrad[,1])+vars$elev*0


####


rr <- rast(vars)
rr$flow_dir <- terrain(rr$elev,"flowdir")
rr$flow_dirladl0 <- flowdirD8ltd(rr$elev,lambda=0,deviation_type="lad") 
##or flowdirD8lad(rr$elev,lambda=0) 

rr$flow_dirladlm <- flowdirD8ltd(rr$elev,lambda=0.5,deviation_type="lad") 
##or flowdirD8lad(rr$elev,lambda=0.5) 

rr$flow_dirladl1 <- flowdirD8ltd(rr$elev,lambda=1,deviation_type="lad")  
##or flowdirD8lad(rr$elev,lambda=1) 

rr$flow_dirltdl0 <- flowdirD8ltd(rr$elev,lambda=0,deviation_type="ltd")
rr$flow_dirltdlm <- flowdirD8ltd(rr$elev,lambda=0.5,deviation_type="ltd")
rr$flow_dirltdl1 <- flowdirD8ltd(rr$elev,lambda=1,deviation_type="ltd")

plot(rr$elev)
arrows_on_rast(rr$flow_dir, unit="flowdir",col="black")

#> NULL
plot(rr$elev)
arrows_on_rast(rr$flow_dirladl0, unit="flowdir",col="blue")

#> NULL
plot(rr$elev)
arrows_on_rast(rr$flow_dirladlm, unit="flowdir",col="black")

#> NULL
plot(rr$elev)
arrows_on_rast(rr$flow_dirladl1, unit="flowdir",col="black")

#> NULL

plot(rr$elev)
arrows_on_rast(rr$flow_dirltdl0, unit="flowdir",col="blue")

#> NULL
plot(rr$elev)
arrows_on_rast(rr$flow_dirltdlm, unit="flowdir",col="black")

#> NULL
plot(rr$elev)
arrows_on_rast(rr$flow_dirltdl1, unit="flowdir",col="black")

#> NULL





### Further examples

elev <- rast(system.file('ex/elev.tif',package="terra"))

flowdirlad <- flowdirD8lad(elev,lambda=0.5)
flowdirlad_copy <- flowdirD8ltd(elev,lambda=0.5,deviation_type="lad")

if (!all(values(flowdirlad_copy==flowdirlad),na.rm=TRUE)) {
  stop("Something in flowdirD8lad  went wrong!") 
}

flowdirladl0 <- flowdirD8lad(elev,lambda=0)
flowdirladl0_copy <- flowdirD8ltd(elev,lambda=0,deviation_type="lad")

if (!all(values(flowdirladl0_copy==flowdirladl0),na.rm=TRUE)) {
  stop("Something in flowdirD8lad  went wrong!") 
}

flowdir <- terrain(elev,"flowdir")
flowdir[flowdirladl0==0] <- 0 
flowdir==flowdirladl0
#> class       : SpatRaster
#> size        : 90, 95, 1  (nrow, ncol, nlyr)
#> resolution  : 0.008333333, 0.008333333  (x, y)
#> extent      : 5.741667, 6.533333, 49.44167, 50.19167  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> source(s)   : memory
#> varname     : elev
#> name        : flowdir
#> min value   :       0
#> max value   :       1



```

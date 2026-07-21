# Path-Based Nondisperive Flow Direction

Computes nondispersive flow/drainage direction according to path-based
methods over grid-based digital elevation models. This is an alternative
to
[`terrain`](https://rspatial.github.io/terra/reference/terrain.md)`(v="flowdir")`.

## Usage

``` r
# S4 method for class 'SpatRaster'
flowDir(x, lambda=0.5, deviation_type=c("ltd","lad"), max_iters=10^6, filename="", ...)
```

## Arguments

- x:

  SpatRaster with elevation data

- lambda:

  Parameter between 0 and 1

- deviation_type:

  Character. Available options are `"ltd"` (the default) and `"lad"`. If
  `"ltd"`, nondispersive flow directions are determined using the least
  transversal deviation (LTD) criterion. If `"lad"`, nondispersive flow
  directions are determined using the least angular deviation (LAD)
  criterion. See Orlandini et al. (2003) for details.

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

The algorithm is an adaptation of the one proposed by Li et al. (2022)
and Orlandini et al. (2003). This function is experimental and under
development: results are to be verified.

## See also

[`terrain`](https://rspatial.github.io/terra/reference/terrain.md),
[`watershed`](https://rspatial.github.io/terra/reference/watershed.md),
[`flowAccumulation`](https://rspatial.github.io/terra/reference/flowAccumulation.md)

## Author

Emanuele Cordano

## References

Orlandini, S., G. Moretti, M. Franchini, B. Aldighieri, and B. Testa
(2003). Path-based methods for the determination of nondispersive
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

## Elevation raster
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
fdir1 <- flowDir(elev1, lambda=1)
fdir2 <- flowDir(elev2, lambda=1)

elev <- rast(system.file('ex/elev.tif',package="terra"))

fdirlad1 <- flowDir(elev, lambda=0.5, deviation_type="lad")
fdirlad2 <- flowDir(elev, lambda=0.5, deviation_type="lad")

###### EXPLICATIVE LONG EXAMPLE 

library(calculus)
#> 
#> Attaching package: ‘calculus’
#> The following object is masked from ‘package:terra’:
#> 
#>     wrap

color_palette <- terrain.colors(50)
### ---------------------------------------------------------
### 1. GRID DEFINITION
### ---------------------------------------------------------
dx <- dy <- 1
xmin <- -10
xmax <- 10
ymin <- -10
ymax <- 10

x <- seq(xmin, xmax, by = dx)
y <- seq(ymin, ymax, by = dy)

vars <- expand.grid(x = x, y = y)

### ---------------------------------------------------------
### 2. DEFINE TWO ANALYTICAL ELEVATION FUNCTIONS
### ---------------------------------------------------------

model=c("vshape","paraboloid","planar")[1]

if (model=="vshape") {
  elev_f1 <-  "2.0*x*tanh(100*x) + 0.50*y"  
  elev_f2 <-   "2.0*x*tanh(100*x) - 0.50*y"
} else if ((model=="paraboloid")) {
  elev_f1 <- " 10*(x^2+y^2)^0.8"
  elev_f2 <- " -10*(x^2+y^2)^0.8"
} else if ((model=="planar")) {
  elev_f1 <- " 2.0*x+0.225*y"
  elev_f2 <- "2.0*x-0.225*y"
}else {
  elev_f1 <- NULL
  elev_f2 <- NULL
}

### ---------------------------------------------------------
### 3. EVALUATE BOTH ELEVATION MODELS
### ---------------------------------------------------------

vars$elev1 <- evaluate(elev_f1, var = vars[,c("x","y")],params=list(ymin=ymin,ymax=ymax))
vars$elev2 <- evaluate(elev_f2, var = vars[,c("x","y")],params=list(ymin=ymin,ymax=ymax))

### ---------------------------------------------------------
### 5. BUILD RASTERS
### ---------------------------------------------------------

r1 <- rast(vars[,c("x","y","elev1")], type="xyz")
r2 <- rast(vars[,c("x","y","elev2")], type="xyz")

names(r1) <- "elev"
names(r2) <- "elev"

### ---------------------------------------------------------
### 6. FLOW DIRECTIONS (terra + LAD/LTD)
### ---------------------------------------------------------
deviation_types=c("lad","ltd")
lambdas <- c(0,0.5,1)

# Elevation Model A
r1$flow_d8 <- terrain(r1$elev, "flowdir")
# Elevation Model B
r2$flow_d8 <- terrain(r2$elev, "flowdir")

for (lambda in lambdas) {
  for (deviation_type in deviation_types) {
    it <- sprintf("flow_%s%03d", deviation_type, round(lambda*10))
    r1[[it]] <- flowDir(r1$elev, lambda=lambda, deviation_type=deviation_type)
    r2[[it]] <- flowDir(r2$elev, lambda=lambda, deviation_type=deviation_type)
  }
}

### ---------------------------------------------------------
### 7. VISUALIZATION EXAMPLES
### ---------------------------------------------------------
nrr <- 1+length(lambdas)*length(deviation_types)
par(mfrow=c(nrr,2))

# Model A
plot(r1$elev, main="Elevation Model A", col=color_palette)
as.arrows(r1$flow_d8, unit="flowdir", col="black")

# Model B
plot(r2$elev, main="Elevation Model B", col=color_palette)
as.arrows(r2$flow_d8, unit="flowdir", col="black")

for (lambda in lambdas) {
  for (deviation_type in deviation_types) {
    ita <- sprintf("flow_%s%03d", deviation_type, round(lambda*10))
    plot(r1$elev, main=sprintf("Model A - %s lambda=%1.2f",toupper(deviation_type),lambda),
         col=color_palette)
    terra::as.arrows(r1[[ita]], unit="flowdir", col="black")
    
    itb <- sprintf("flow_%s%03d", deviation_type, round(lambda*10))
    plot(r2$elev, main=sprintf("Model B - %s lambda=%1.2f",toupper(deviation_type),lambda),
         col=color_palette)
    terra::as.arrows(r2[[itb]], unit="flowdir", col="black")
  }
}


par(mfrow=c(1,1))

message("Analysis complete")
#> Analysis complete





```

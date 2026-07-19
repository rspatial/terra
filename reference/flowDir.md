# Path-Based Nondisperive Flow Direction

Computes non-disperive flow/draiage direction according to path-based
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

  Character. Default is the first element of c("ltd", "lad"). If `"ltd"`
  (the default), nondispersive flow directions are determined using the
  least transversal deviation criterion (LTD). If `"lad"`, nondispersive
  flow directions are determined using the least angular deviation (LAD)
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

The algorithm is an adaptation of the one proposed by Li et al, 2021 and
Orlandini et al, 2003 This function is experimental and under
development: results are to be verified.

## See also

[`terrain`](https://rspatial.github.io/terra/reference/terrain.md),
[`watershed`](https://rspatial.github.io/terra/reference/watershed.md),
[`flowAccumulation`](https://rspatial.github.io/terra/reference/flowAccumulation.md)

## Author

Emanuele Cordano

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
fdir2 <- flowDir(elev2,lambda=1)

elev <- rast(system.file('ex/elev.tif',package="terra"))

fdirlad1 <- flowDir(elev, lambda=0.5, deviation_type="lad")
fdirlad2 <- flowDir(elev, lambda=0.5, deviation_type="lad")
```

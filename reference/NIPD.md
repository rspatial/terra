# Number of immediate adjacent cells flowing into each cell

Compute the number of immediate adjacent cells flowing into each cell

## Usage

``` r
# S4 method for class 'SpatRaster'
NIDP(x, filename="",...)
```

## Arguments

- x:

  SpatRaster with flow-direction. see
  [`terrain`](https://rspatial.github.io/terra/reference/terrain.md)

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## Details

NDIP is computed first to compute flow-accumulation with the algorithm
by Zhou at al, 2019.

## References

Zhou, G., Wei, H. & Fu, S. A fast and simple algorithm for calculating
flow accumulation matrices from raster digital elevation. Front. Earth
Sci. 13, 317–326 (2019). https://doi.org/10.1007/s11707-018-0725-9
<https://link.springer.com/article/10.1007/s11707-018-0725-9>

## See also

[`flowAccumulation`](https://rspatial.github.io/terra/reference/flowAccumulation.md)

## Author

Emanuele Cordano

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
    elev2[r,c] <- 10+5*(abs(x))-0.001*y ### 5*(x^2+y^2)
  }
} 


## Elevation Raster 
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


## Flow Direction Raster
flowdir1<- terrain(elev1,v="flowdir")
flowdir2<- terrain(elev2,v="flowdir")


t(array(flowdir1[],rev(dim(flowdir1)[1:2])))
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    2    2    2    4    4    4    8    8    8
#>  [2,]    2    2    2    4    4    4    8    8    8
#>  [3,]    2    2    2    4    4    4    8    8    8
#>  [4,]    1    1    1    2    4    8   16   16   16
#>  [5,]    1    1    1    1    0   16   16   16   16
#>  [6,]    1    1    1  128   64   32   16   16   16
#>  [7,]  128  128  128   64   64   64   32   32   32
#>  [8,]  128  128  128   64   64   64   32   32   32
#>  [9,]  128  128  128   64   64   64   32   32   32
t(array(flowdir2[],rev(dim(flowdir2)[1:2])))
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

plot(flowdir1)

plot(flowdir2)


## 
nidp1 <- NIDP((flowdir1))
nidp2 <- NIDP((flowdir2))

t(array(nidp1[],rev(dim(nidp1)[1:2])))
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    0    0    0    0    0    0    0    0    0
#>  [2,]    0    1    1    2    1    2    1    1    0
#>  [3,]    0    1    1    2    1    2    1    1    0
#>  [4,]    0    2    2    3    1    3    2    2    0
#>  [5,]    0    1    1    1    9    1    1    1    0
#>  [6,]    0    2    2    3    1    3    2    2    0
#>  [7,]    0    1    1    2    1    2    1    1    0
#>  [8,]    0    1    1    2    1    2    1    1    0
#>  [9,]    0    0    0    0    0    0    0    0    0
t(array(nidp2[],rev(dim(nidp2)[1:2])))
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    0    1    1    1    2    1    1    1    0
#>  [2,]    0    1    1    1    3    1    1    1    0
#>  [3,]    0    1    1    1    3    1    1    1    0
#>  [4,]    0    1    1    1    3    1    1    1    0
#>  [5,]    0    1    1    1    3    1    1    1    0
#>  [6,]    0    1    1    1    3    1    1    1    0
#>  [7,]    0    1    1    1    3    1    1    1    0
#>  [8,]    0    1    1    1    3    1    1    1    0
#>  [9,]    0    1    1    1    4    1    1    1    0

plot(nidp1)

plot(nidp2)

```

# Get the coordinates of SpatVector geometries or SpatRaster cells

Get the coordinates of a SpatVector or SpatRaster cells. A matrix or
data.frame of the x (longitude) and y (latitude) coordinates is
returned.

## Usage

``` r
# S4 method for class 'SpatVector'
crds(x, df=FALSE, list=FALSE)

# S4 method for class 'SpatRaster'
crds(x, df=FALSE, na.rm=TRUE, na.all=FALSE)
```

## Arguments

- x:

  SpatRaster or SpatVector

- df:

  logical. If `TRUE` a `data.frame` is returned instead of a matrix

- list:

  logical. If `TRUE` a `list` is returned instead of a matrix

- na.rm:

  logical. If `TRUE` cells that are `NA` are excluded. Ignored if the
  SpatRaster is a template with no associated cell values

- na.all:

  logical. If `TRUE` cells are only ignored if `na.rm=TRUE` and their
  value is `NA` for **all** layers instead of for `any` layer

## Value

matrix or data.frame

## See also

[`geom`](https://rspatial.github.io/terra/reference/geometry.md) returns
the complete structure of SpatVector geometries. For SpatRaster see
[`xyFromCell`](https://rspatial.github.io/terra/reference/xyCellFrom.md)

## Examples

``` r
x1 <- rbind(c(-175,-20), c(-140,55), c(10, 0), c(-140,-60))
x2 <- rbind(c(-125,0), c(0,60), c(40,5), c(15,-45))
x3 <- rbind(c(-10,0), c(140,60), c(160,0), c(140,-55))
x4 <- rbind(c(80,0), c(105,13), c(120,2), c(105,-13))
z <- rbind(cbind(object=1, part=1, x1), cbind(object=2, part=1, x2), 
           cbind(object=3, part=1, x3), cbind(object=3, part=2, x4))
colnames(z)[3:4] <- c('x', 'y')
z <- cbind(z, hole=0)
z[(z[, "object"]==3 & z[,"part"]==2), "hole"] <- 1

p <- vect(z, "polygons")
crds(p)
#>          x   y
#>  [1,] -175 -20
#>  [2,] -140  55
#>  [3,]   10   0
#>  [4,] -140 -60
#>  [5,] -175 -20
#>  [6,] -125   0
#>  [7,]    0  60
#>  [8,]   40   5
#>  [9,]   15 -45
#> [10,] -125   0
#> [11,]  -10   0
#> [12,]  140  60
#> [13,]  160   0
#> [14,]  140 -55
#> [15,]  -10   0
#> [16,]   80   0
#> [17,]  105  13
#> [18,]  120   2
#> [19,]  105 -13
#> [20,]   80   0

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
g <- crds(v)
head(g)
#>             x        y
#> [1,] 6.026519 50.17767
#> [2,] 6.031361 50.16563
#> [3,] 6.035646 50.16410
#> [4,] 6.042747 50.16157
#> [5,] 6.043894 50.16116
#> [6,] 6.048243 50.16008
```

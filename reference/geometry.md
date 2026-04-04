# Get the geometry (coordinates) of a SpatVector

Get the geometry of a SpatVector. If `wkt=FALSE`, this is a five-column
matrix or data.frame: the vector object ID, the IDs for the parts of
each object (e.g. five polygons that together are one spatial object),
the x (longitude) and y (latitude) coordinates, and a flag indicating
whether the part is a "hole" (only relevant for polygons).

If `wkt=TRUE`, the "well-known text" representation is returned as a
character vector. If `hex=TRUE`, the "hexadecimal" representation is
returned as a character vector. If `wkb=TRUE`, the "well-known binary"
representation is returned as a list of raw vectors.

## Usage

``` r
# S4 method for class 'SpatVector'
geom(x, wkt=FALSE, hex=FALSE, wkb=FALSE, df=FALSE, list=FALSE, xnm="x", ynm="y")
```

## Arguments

- x:

  SpatVector

- wkt:

  logical. If `TRUE` the WKT geometry is returned (unless `hex` is also
  `TRUE`)

- hex:

  logical. If `TRUE` the hexadecimal geometry is returned

- wkb:

  logical. If `TRUE` the raw WKB geometry is returned (unless either of
  `hex` or `wkt` is also `TRUE`)

- df:

  logical. If `TRUE` a `data.frame` is returned instead of a matrix
  (only if `wkt=FALSE`, `hex=FALSE`, and `list=FALSE`)

- list:

  logical. If `TRUE` a nested `list` is returned with data.frames of
  coordinates

- xnm:

  character. If `list=TRUE` the "x" column name for the coordinates
  data.frame

- ynm:

  character. If `list=TRUE` the "y" column name for the coordinates
  data.frame

## Value

matrix, vector, data.frame, or list

## See also

[`crds`](https://rspatial.github.io/terra/reference/crds.md),
[`xyFromCell`](https://rspatial.github.io/terra/reference/xyCellFrom.md)

## Examples

``` r
x1 <- rbind(c(-175,-20), c(-140,55), c(10, 0), c(-140,-60))
x2 <- rbind(c(-125,0), c(0,60), c(40,5), c(15,-45))
x3 <- rbind(c(-10,0), c(140,60), c(160,0), c(140,-55))
x4 <- rbind(c(80,0), c(105,13), c(120,2), c(105,-13))
z <- rbind(cbind(object=1, part=1, x1), cbind(object=2, part=1, x2), 
           cbind(object=3, part=1, x3), cbind(object=3, part=2,  x4))
colnames(z)[3:4] <- c('x', 'y')
z <- cbind(z, hole=0)
z[(z[, "object"]==3 & z[,"part"]==2), "hole"] <- 1

p <- vect(z, "polygons")
geom(p)
#>       geom part    x   y hole
#>  [1,]    1    1 -175 -20    0
#>  [2,]    1    1 -140  55    0
#>  [3,]    1    1   10   0    0
#>  [4,]    1    1 -140 -60    0
#>  [5,]    1    1 -175 -20    0
#>  [6,]    2    1 -125   0    0
#>  [7,]    2    1    0  60    0
#>  [8,]    2    1   40   5    0
#>  [9,]    2    1   15 -45    0
#> [10,]    2    1 -125   0    0
#> [11,]    3    1  -10   0    0
#> [12,]    3    1  140  60    0
#> [13,]    3    1  160   0    0
#> [14,]    3    1  140 -55    0
#> [15,]    3    1  -10   0    0
#> [16,]    3    1   80   0    1
#> [17,]    3    1  105  13    1
#> [18,]    3    1  120   2    1
#> [19,]    3    1  105 -13    1
#> [20,]    3    1   80   0    1

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
g <- geom(v)
head(g)
#>      geom part        x        y hole
#> [1,]    1    1 6.026519 50.17767    0
#> [2,]    1    1 6.031361 50.16563    0
#> [3,]    1    1 6.035646 50.16410    0
#> [4,]    1    1 6.042747 50.16157    0
#> [5,]    1    1 6.043894 50.16116    0
#> [6,]    1    1 6.048243 50.16008    0

w <- geom(v, wkt=TRUE)
substr(w, 1, 60)
#>  [1] "POLYGON ((6.0265188199999997 50.1776695300000029, 6.03136109"
#>  [2] "POLYGON ((6.1783676099999996 49.8768234300000017, 6.18547869"
#>  [3] "POLYGON ((5.8813781699999996 49.8701477099999977, 5.88167238"
#>  [4] "POLYGON ((6.1313085599999999 49.9725646999999995, 6.13429117"
#>  [5] "POLYGON ((5.9779286400000000 50.0260162399999970, 5.98231220"
#>  [6] "POLYGON ((6.3855319000000001 49.8370285000000024, 6.38859987"
#>  [7] "POLYGON ((6.3166646999999996 49.6233749399999979, 6.31834983"
#>  [8] "POLYGON ((6.4251575499999998 49.7316436799999977, 6.42656994"
#>  [9] "POLYGON ((5.9983120000000003 49.6999244699999991, 5.99863242"
#> [10] "POLYGON ((6.0394744899999999 49.4482612600000024, 6.03690576"
#> [11] "POLYGON ((6.1559634200000000 49.6850471500000026, 6.15928363"
#> [12] "POLYGON ((6.0679822000000003 49.8284645100000034, 6.07192183"
```

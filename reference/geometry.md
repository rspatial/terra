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
#>  [1] "POLYGON ((6.026519 50.17767, 6.031361 50.165627, 6.035646 50"
#>  [2] "POLYGON ((6.178368 49.876823, 6.185479 49.870525, 6.189417 4"
#>  [3] "POLYGON ((5.881378 49.870148, 5.881672 49.868866, 5.886637 4"
#>  [4] "POLYGON ((6.131309 49.972565, 6.134291 49.972382, 6.139316 4"
#>  [5] "POLYGON ((5.977929 50.026016, 5.982312 50.022949, 5.981743 5"
#>  [6] "POLYGON ((6.385532 49.837029, 6.3886 49.833683, 6.390184 49."
#>  [7] "POLYGON ((6.316665 49.623375, 6.31835 49.623157, 6.320131 49"
#>  [8] "POLYGON ((6.425158 49.731644, 6.42657 49.73082, 6.427332 49."
#>  [9] "POLYGON ((5.998312 49.699924, 5.998632 49.698559, 5.998956 4"
#> [10] "POLYGON ((6.039474 49.448261, 6.036906 49.448696, 6.036822 4"
#> [11] "POLYGON ((6.155963 49.685047, 6.159284 49.685036, 6.161457 4"
#> [12] "POLYGON ((6.067982 49.828465, 6.071922 49.825478, 6.073236 4"
```

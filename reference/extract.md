# Extract values from a SpatRaster

Extract values from a SpatRaster for a set of locations. The locations
can be a SpatVector (points, lines, polygons), a data.frame or matrix
with (x, y) or (longitude, latitude – in that order!) coordinates, or a
vector with cell numbers.

When argument `y` is a `SpatVector` the first column has the ID (record
number) of the `SpatVector` used (unless you set `ID=FALSE`).

Alternatively, you can use
[`zonal`](https://rspatial.github.io/terra/reference/zonal.md) after
using
[`rasterize`](https://rspatial.github.io/terra/reference/rasterize.md)
with a `SpatVector` (this may be more efficient in some cases).

## Usage

``` r
# S4 method for class 'SpatRaster,SpatVector'
extract(x, y, fun=NULL, method="simple", cells=FALSE, xy=FALSE,
    ID=TRUE, weights=FALSE, exact=FALSE, touches=is.lines(y), small=TRUE,
  layer=NULL, bind=FALSE, raw=FALSE, search_radius=0, ...)

# S4 method for class 'SpatRaster,data.frame'
extract(x, y, ...)

# S4 method for class 'SpatRaster,SpatExtent'
extract(x, y, ...)

# S4 method for class 'SpatRaster,matrix'
extract(x, y, ...)

# S4 method for class 'SpatRaster,numeric'
extract(x, y, ...)


# S4 method for class 'SpatVector,SpatVector'
extract(x, y, count=FALSE)
```

## Arguments

- x:

  SpatRaster or SpatVector of polygons

- y:

  SpatVector (points, lines, or polygons). Alternatively, for points, a
  2-column matrix or data.frame (x, y) or (lon, lat). Or a vector with
  cell numbers

- fun:

  function to summarize the extracted data by line or polygon geometry.
  You can use `fun=table` to tabulate raster values for each line or
  polygon geometry. If `weights=TRUE` or `exact=TRUE` only `mean`,
  `sum`, `min`, `max` and `table` are accepted — and these functions
  will consider the fraction of a cell that is covered when computing
  the mean or the sum). Ignored if `y` has point geometry

- method:

  character. method for extracting values with points ("simple" or
  "bilinear"). With "simple" values for the cell a point falls in are
  returned. With "bilinear" the returned values are interpolated from
  the values of the four nearest raster cells

- cells:

  logical. If `TRUE` the cell numbers are also returned, unless `fun` is
  not `NULL`. Also see
  [`cells`](https://rspatial.github.io/terra/reference/cells.md)

- xy:

  logical. If `TRUE` the coordinates of the cells are also returned,
  unless `fun` is not `NULL`. See
  [`xyFromCell`](https://rspatial.github.io/terra/reference/xyCellFrom.md)

- ID:

  logical. Should an ID column be added? If so, the first column
  returned has the IDs (record numbers) of `y`

- weights:

  logical. If `TRUE` and `y` has polygons, the approximate fraction of
  each cell that is covered is returned as well. This changes the effect
  of argument `fun`

- exact:

  logical. If `TRUE` and `y` has polygons, the exact fraction of each
  cell that is covered is returned as well. This changes the effect of
  argument `fun`

- touches:

  logical. If `TRUE`, values for all cells touched by lines or polygons
  are extracted, not just those on the line render path, or whose center
  point is within the polygon. Not relevant for points; and always
  considered `TRUE` when `weights=TRUE` or `exact=TRUE`

- small:

  logical. If `TRUE`, values for all cells in touched polygons are
  extracted if none of the cells center points is within the polygon;
  even if `touches=FALSE`

- layer:

  character or numeric to select the layer to extract from for each
  geometry. If `layer` is a character it can be a name in `y` or a
  vector of layer names. If it is numeric, it must be integer values
  between `1` and `nlyr(x)`

- bind:

  logical. If `TRUE`, a SpatVector is returned consisting of the input
  SpatVector `y` and the `cbind`-ed extracted values

- raw:

  logical. If `TRUE`, a matrix is returned with the "raw" numeric cell
  values. If `FALSE`, a data.frame is returned and the cell values are
  transformed to factor, logical, or integer values, where appropriate

- search_radius:

  positive number. A search-radius that is used when `y` has point
  geometry. If this value is larger than zero, it is the maximum
  distance used to find the a cell with a value that is nearest to the
  cell that the point falls in if that cell that has a missing (`NA`)
  value. The value of this nearest cell, the distance to the original
  cell, and the new cell number are returned. The radius should be
  expressed in m if the data have lon/lat coordinates or in the distance
  unit of the crs in other cases (typically also m). For lon/lat data,
  the mean latitude of the points is used to compute the distances, so
  this may be imprecise for data with a large latitudinal range

- ...:

  additional arguments to `fun` if `y` is a SpatVector. For example
  `na.rm=TRUE`. Or arguments passed to the `SpatRaster,SpatVector`

- count:

  logical. If `TRUE` and `x` has polygons geometry and `y` has points
  geometry, the number of points in polygons is returned

## Value

data.frame, matrix or SpatVector

## See also

[`values`](https://rspatial.github.io/terra/reference/values.md)`, `[`zonal`](https://rspatial.github.io/terra/reference/zonal.md)`, `[`extractAlong`](https://rspatial.github.io/terra/reference/extractAlong.md)`, `[`extractRange`](https://rspatial.github.io/terra/reference/extractRange.md)`, `[`rapp`](https://rspatial.github.io/terra/reference/rapp.md)

## Examples

``` r
r <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
values(r) <- 1:25
xy <- cbind(lon=c(0.5,2.5), lat=c(0.5,2.5))
p <- vect(xy, crs="+proj=longlat +datum=WGS84")

extract(r, xy)
#> Error in extract(r, xy): object 'ID' not found
extract(r, p)
#>   ID lyr.1
#> 1  1    21
#> 2  2    13

r[1,]
#>   lyr.1
#> 1     1
#> 2     2
#> 3     3
#> 4     4
#> 5     5
r[5]
#>   lyr.1
#> 1     5
r[,5]
#>   lyr.1
#> 1     5
#> 2    10
#> 3    15
#> 4    20
#> 5    25

r[c(0:2, 99:101)]
#>   lyr.1
#> 1     1
#> 2     2
#> 3    NA
#> 4    NA
#> 5    NA

f <- system.file("ex/meuse.tif", package="terra")
r <- rast(f)

xy <- cbind(179000, 330000)
xy <- rbind(xy-100, xy, xy+1000)
extract(r, xy)
#> Error in extract(r, xy): object 'ID' not found

p <- vect(xy)
g <- geom(p)
g
#>      geom part      x      y hole
#> [1,]    1    1 178900 329900    0
#> [2,]    2    1 179000 330000    0
#> [3,]    3    1 180000 331000    0

extract(r, p)
#>   ID meuse
#> 1  1   378
#> 2  2   251
#> 3  3   208

x <- r + 10
extract(x, p)
#>   ID meuse
#> 1  1   388
#> 2  2   261
#> 3  3   218

i <- cellFromXY(r, xy)
x[i]
#>   meuse
#> 1   388
#> 2   261
#> 3   218
r[i]
#>   meuse
#> 1   378
#> 2   251
#> 3   208

y <- c(x,x*2,x*3)
y[i]
#>   meuse meuse meuse
#> 1   388   776  1164
#> 2   261   522   783
#> 3   218   436   654

## extract with a polygon
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
v <- v[1:2,]

rf <- system.file("ex/elev.tif", package="terra")
x <- rast(rf)
extract(x, v, mean, na.rm=TRUE)
#>   ID elevation
#> 1  1  467.1052
#> 2  2  333.8629

z <- rast(v, resolution=.1, names="test")
values(z) <- 1:ncell(z)
e <- extract(z, v, ID=TRUE)
e
#>   ID test
#> 1  1    2
#> 2  1    3
#> 3  1    6
#> 4  1    7
#> 5  1    8
#> 6  2   13
#> 7  2   17
#> 8  2   18
#> 9  2   19
tapply(e[,2], e[,1], mean, na.rm=TRUE)
#>     1     2 
#>  5.20 16.75 

x <- c(z, z*2, z/3)
names(x) <- letters[1:3]

e <- extract(x, v, ID=TRUE)
de <- data.frame(e)
aggregate(de[,2:4], de[,1,drop=FALSE], mean)
#>   ID     a    b        c
#> 1  1  5.20 10.4 1.733333
#> 2  2 16.75 33.5 5.583333
```

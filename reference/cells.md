# Get cell numbers

Get the cell numbers covered by a SpatVector or SpatExtent. Or that
match values in a vector; or all non `NA` values.

## Usage

``` r
# S4 method for class 'SpatRaster,missing'
cells(x, y)

# S4 method for class 'SpatRaster,numeric'
cells(x, y, pairs=FALSE)

# S4 method for class 'SpatRaster,SpatVector'
cells(x, y, method="simple", weights=FALSE, exact=FALSE, 
        touches=is.lines(y), small=TRUE)

# S4 method for class 'SpatRaster,SpatExtent'
cells(x, y)
```

## Arguments

- x:

  SpatRaster

- y:

  SpatVector, SpatExtent, 2-column matrix representing points, numeric
  representing values to match, or missing

- method:

  character. Method for getting cell numbers for points. The default is
  "simple", the alternative is "bilinear". If it is "bilinear", the four
  nearest cells and their weights are returned

- weights:

  logical. If `TRUE` and `y` has polygons, the approximate fraction of
  each cell that is covered is returned as well

- pairs:

  logical. If `TRUE` the cell values matched area also returned

- exact:

  logical. If `TRUE` and `y` has polygons, the exact fraction of each
  cell that is covered is returned as well

- touches:

  logical. If `TRUE`, values for all cells touched by lines or polygons
  are extracted, not just those on the line render path, or whose center
  point is within the polygon. Not relevant for points

- small:

  logical. If `TRUE`, values for all cells in touched polygons are
  extracted if none of the cells center points is within the polygon;
  even if `touches=FALSE`

## Value

numeric vector or matrix

## Examples

``` r
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
r[c(1:25, 31:100)] <- NA
r <- ifel(r > 28, r + 10, r)

# all cell numbers of cells that are not NA
cells(r)
#> [1] 26 27 28 29 30

# cell numbers that match values
x <- cells(r, c(28,38))
x$lyr.1
#> [1] 28

# cells for points
m <- cbind(x=c(0,10,-30), y=c(40,-10,20))
cellFromXY(r, m)
#> [1] 26 56 35

v <- vect(m)
cells(r, v)
#>      ID cell
#> [1,]  1   26
#> [2,]  2   56
#> [3,]  3   35
cells(r, v, method="bilinear")
#>      ID c1 c2 c3 c4        w1        w2         w3         w4
#> [1,]  1 25 26 35 36 0.3611111 0.3611111 0.13888889 0.13888889
#> [2,]  2 55 56 65 66 0.2098765 0.7345679 0.01234568 0.04320988
#> [3,]  3 34 35 44 45 0.2037037 0.4074074 0.12962963 0.25925926

# cells for polygons
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
r <- rast(v)
cv <- cells(r, v) 
```

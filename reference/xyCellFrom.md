# Coordinates from a row, column or cell number and vice versa

Get coordinates of the center of raster cells for a row, column, or cell
number of a SpatRaster. Or get row, column, or cell numbers from
coordinates or from each other.

Cell numbers start at 1 in the upper left corner, and increase from left
to right, and then from top to bottom. The last cell number equals the
number of cells of the SpatRaster (see
[`ncell`](https://rspatial.github.io/terra/reference/dimensions.md)).
Row numbers start at 1 at the top, column numbers start at 1 at the
left.

When computing row, column, or cell numbers from coordinates, and
coordinates fall on the edge of two or four cells, they are assigned to
the right-most and/or lowest cell. That is, in these cases of ambiguity,
the highest row, column, or cell number is returned.

## Usage

``` r
# S4 method for class 'SpatRaster,numeric'
xFromCol(object, col)

# S4 method for class 'SpatRaster,numeric'
yFromRow(object, row)

# S4 method for class 'SpatRaster,numeric'
xyFromCell(object, cell)

# S4 method for class 'SpatRaster,numeric'
xFromCell(object, cell)

# S4 method for class 'SpatRaster,numeric'
yFromCell(object, cell)

# S4 method for class 'SpatRaster,numeric'
colFromX(object, x)

# S4 method for class 'SpatRaster,numeric'
rowFromY(object, y)

# S4 method for class 'SpatRaster,numeric,numeric'
cellFromRowCol(object, row, col)

# S4 method for class 'SpatRaster,numeric,numeric'
cellFromRowColCombine(object, row, col)

# S4 method for class 'SpatRaster,numeric,numeric'
rowColCombine(object, row, col)

# S4 method for class 'SpatRaster,numeric'
rowFromCell(object, cell)

# S4 method for class 'SpatRaster,numeric'
colFromCell(object, cell)

# S4 method for class 'SpatRaster,numeric'
rowColFromCell(object, cell)

# S4 method for class 'SpatRaster,matrix'
cellFromXY(object, xy)
```

## Arguments

- object:

  SpatRaster

- cell:

  integer. cell number(s)

- col:

  integer. column number(s) or missing (equivalent to all columns)

- row:

  integer. row number(s) or missing (equivalent to all rows)

- x:

  x coordinate(s)

- y:

  y coordinate(s)

- xy:

  matrix of x and y coordinates

## Value

xFromCol, yFromCol, xFromCell, yFromCell: vector of x or y coordinates

xyFromCell: matrix(x,y) with coordinate pairs

colFromX, rowFromY, cellFromXY, cellFromRowCol, rowFromCell,
colFromCell: vector of row, column, or cell numbers

rowColFromCell, rowColCombine: matrix of row and column numbers

## See also

[`crds`](https://rspatial.github.io/terra/reference/crds.md)

## Examples

``` r
r <- rast()

xFromCol(r, c(1, 120, 180))
#> [1] -179.5  -60.5   -0.5
yFromRow(r, 90)
#> [1] 0.5
xyFromCell(r, 10000)
#>         x    y
#> [1,] 99.5 62.5
xyFromCell(r, c(0, 1, 32581, ncell(r), ncell(r)+1))
#>           x     y
#> [1,]    NaN   NaN
#> [2,] -179.5  89.5
#> [3,]    0.5  -0.5
#> [4,]  179.5 -89.5
#> [5,]    NaN   NaN

cellFromRowCol(r, 5, 5)
#> [1] 1445
cellFromRowCol(r, 1:2, 1:2)
#> [1]   1 362
cellFromRowCol(r, 1, 1:3)
#> [1] 1 2 3

# all combinations
cellFromRowColCombine(r, 1:2, 1:2)
#> [1]   1   2 361 362

colFromX(r, 10)
#> [1] 191
rowFromY(r, 10)
#> [1] 81
xy <- cbind(lon=c(10,5), lat=c(15, 88))
cellFromXY(r, xy)
#> [1] 27191   906

# if no row/col specified all are returned
range(xFromCol(r))
#> [1] -179.5  179.5
length(yFromRow(r))
#> [1] 180
```

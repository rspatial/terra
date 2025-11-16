# Adjacent cells or polygons

Identify cells that are adjacent to a set of raster cells. Or identify
adjacent polygons

## Usage

``` r
# S4 method for class 'SpatRaster'
adjacent(x, cells, directions="rook", pairs=FALSE, include=FALSE, symmetrical=FALSE)

# S4 method for class 'SpatVector'
adjacent(x, type="rook", pairs=TRUE, symmetrical=FALSE)
```

## Arguments

- x:

  SpatRaster, or SpatVector of polygons

- cells:

  vector of cell numbers for which adjacent cells should be found. Cell
  numbers start with 1 in the upper-left corner and increase from left
  to right and from top to bottom

- directions:

  character or matrix to indicated the directions in which cells are
  considered connected. The following character values are allowed:
  "rook" or "4" for the horizontal and vertical neighbors; "bishop" to
  get the diagonal neighbors; "queen" or "8" to get the vertical,
  horizontal and diagonal neighbors; or "16" for knight and one-cell
  queen move neighbors. If `directions` is a matrix it should have odd
  dimensions and have logical (or 0, 1) values

- pairs:

  logical. If `TRUE`, a two-column matrix of pairs of adjacent cells is
  returned. If `x` is a `SpatRaster` and `pairs` is `FALSE`, an `n*m`
  matrix is returned where the number of rows `n` is `length(cells)` and
  the number of columns `m` is the number of neighbors requested with
  `directions`

- include:

  logical. Should the focal cells be included in the result?

- type:

  character. One of "rook", "queen", "touches", or "intersects". "queen"
  and "touches" are synonyms. "rook" exclude polygons that touch at a
  single node only. "intersects" includes polygons that touch or overlap

- symmetrical:

  logical. If `TRUE` and `pairs=TRUE`, an adjacent pair is only included
  once. For example, if polygon 1 is adjacent to polygon 3, the implied
  adjacency between 3 and 1 is not reported

## Note

When using global lon/lat rasters, adjacent cells at the other side of
the date-line are included.

## See also

[`relate`](https://rspatial.github.io/terra/reference/relate.md)`, `[`nearby`](https://rspatial.github.io/terra/reference/nearby.md)`, `[`nearest`](https://rspatial.github.io/terra/reference/nearby.md)

## Value

matrix

## Examples

``` r
r <- rast(nrows=10, ncols=10)
adjacent(r, cells=c(1, 5, 55), directions="queen") 
#>    [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> 0   NaN  NaN  NaN   10    2   20   11   12
#> 4   NaN  NaN  NaN    4    6   14   15   16
#> 54   44   45   46   54   56   64   65   66
r <- rast(nrows=10, ncols=10, crs="+proj=utm +zone=1 +datum=WGS84")
adjacent(r, cells=11, directions="rook") 
#>    [,1] [,2] [,3] [,4]
#> 10    1  NaN   12   21

#same as 
rk <- matrix(c(0,1,0,1,0,1,0,1,0), 3, 3)
adjacent(r, cells=11, directions=rk) 
#>    [,1] [,2] [,3] [,4]
#> 10    1  NaN   12   21

## note that with global lat/lon data the E and W connect
r <- rast(nrows=10, ncols=10, crs="+proj=longlat +datum=WGS84")
adjacent(r, cells=11, directions="rook") 
#>    [,1] [,2] [,3] [,4]
#> 10    1   20   12   21

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
a <- adjacent(v, symmetrical=TRUE)
head(a)
#>      from to
#> [1,]    1  2
#> [2,]    1  4
#> [3,]    1  5
#> [4,]    2  3
#> [5,]    2  4
#> [6,]    2  5
```

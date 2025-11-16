# Coercion to vector, matrix or array

Coercion of a SpatRaster to a vector, matrix or array. Or coerce a
SpatExtent to a vector or matrix

## Usage

``` r
# S4 method for class 'SpatRaster'
as.vector(x, mode='any')

# S4 method for class 'SpatRaster'
as.matrix(x, wide=FALSE, ...)

# S4 method for class 'SpatRaster'
as.array(x)

# S4 method for class 'SpatRasterDataset'
as.array(x)

# S4 method for class 'SpatExtent'
as.vector(x, mode='any')

# S4 method for class 'SpatExtent'
as.matrix(x, ...)
```

## Arguments

- x:

  SpatRaster or SpatVector

- wide:

  logical. If `FALSE` each layer in the SpatRaster becomes a column in
  the matrix and each cell in the SpatRaster becomes a row. If `TRUE`
  each row in the SpatRaster becomes a row in the matrix and each column
  in the SpatRaster becomes a column in the matrix

- mode:

  this argument is ignored

- ...:

  additional arguments (none implemented)

## Value

vector, matrix, or array

## See also

[`as.data.frame`](https://rspatial.github.io/terra/reference/as.data.frame.md)
and
[`as.polygons`](https://rspatial.github.io/terra/reference/as.polygons.md)

## Examples

``` r
r <- rast(ncols=2, nrows=2)
values(r) <- 1:ncell(r)

as.vector(r)
#> [1] 1 2 3 4
as.matrix(r)
#>      lyr.1
#> [1,]     1
#> [2,]     2
#> [3,]     3
#> [4,]     4
as.matrix(r, wide=TRUE)
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    3    4
as.data.frame(r, xy=TRUE)
#>     x   y lyr.1
#> 1 -90  45     1
#> 2  90  45     2
#> 3 -90 -45     3
#> 4  90 -45     4
as.array(r)
#> , , 1
#> 
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    3    4
#> 

as.vector(ext(r))
#> xmin xmax ymin ymax 
#> -180  180  -90   90 
as.matrix(ext(r))
#>       min max
#> [1,] -180 180
#> [2,]  -90  90
```

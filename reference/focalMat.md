# Focal weights matrix

Make a focal ("moving window") weight matrix for use in the
[`focal`](https://rspatial.github.io/terra/reference/focal.md) function.
The sum of the values adds up to one.

## Usage

``` r
focalMat(x, d, type=c('circle', 'Gauss', 'rectangle'), fillNA=FALSE)
```

## Arguments

- x:

  SpatRaster

- d:

  numeric. If `type=circle`, the radius of the circle (in units of the
  crs). If `type=rectangle` the dimension of the rectangle (one or two
  numbers). If `type=Gauss` the size of sigma, and optionally another
  number to determine the size of the matrix returned (default is
  3\*sigma)

- type:

  character indicating the type of filter to be returned

- fillNA:

  logical. If `TRUE`, zeros are set to `NA` such that they are ignored
  in the computations. Only applies to `type="circle"`

## Value

matrix that can be used with
[`focal`](https://rspatial.github.io/terra/reference/focal.md)

## Examples

``` r
r <- rast(ncols=180, nrows=180, xmin=0)
focalMat(r, 2, "circle")
#>            [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,] 0.00000000 0.00000000 0.07692308 0.00000000 0.00000000
#> [2,] 0.00000000 0.07692308 0.07692308 0.07692308 0.00000000
#> [3,] 0.07692308 0.07692308 0.07692308 0.07692308 0.07692308
#> [4,] 0.00000000 0.07692308 0.07692308 0.07692308 0.00000000
#> [5,] 0.00000000 0.00000000 0.07692308 0.00000000 0.00000000

focalMat(r, c(2,3), "rect")
#>            [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,] 0.02857143 0.02857143 0.02857143 0.02857143 0.02857143
#> [2,] 0.02857143 0.02857143 0.02857143 0.02857143 0.02857143
#> [3,] 0.02857143 0.02857143 0.02857143 0.02857143 0.02857143
#> [4,] 0.02857143 0.02857143 0.02857143 0.02857143 0.02857143
#> [5,] 0.02857143 0.02857143 0.02857143 0.02857143 0.02857143
#> [6,] 0.02857143 0.02857143 0.02857143 0.02857143 0.02857143
#> [7,] 0.02857143 0.02857143 0.02857143 0.02857143 0.02857143

# Gaussian filter for square cells
gf <- focalMat(r, 1, "Gauss")
```

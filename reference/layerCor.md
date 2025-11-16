# Correlation and (weighted) covariance

Compute correlation, (weighted) covariance, or similar summary
statistics that compare the values of all pairs of the layers of a
SpatRaster.

## Usage

``` r
# S4 method for class 'SpatRaster'
layerCor(x, fun, w, asSample=TRUE, use="everything", maxcell=Inf, ...)
```

## Arguments

- x:

  SpatRaster

- fun:

  character. The statistic to compute: either "cov" (covariance),
  "weighted.cov" (weighted covariance), or "cor" (pearson correlation
  coefficient). You can also supply your own function that takes two
  vectors as argument to compute a single number

- w:

  SpatRaster with the weights to compute the weighted covariance. It
  should have a single layer and the same geometry as `x`

- asSample:

  logical. If `TRUE`, the statistic for a sample (denominator is `n-1`)
  is computed, rather than for the population (denominator is `n`). Only
  for the standard functions

- use:

  character. To decide how to handle missing values. This must be (an
  abbreviation of) one of "everything", "complete.obs",
  "pairwise.complete.obs", "masked.complete". With
  "pairwise.complete.obs", the value for a pair of layers is computed
  for all cells that are not `NA` in that pair. Therefore, it may be
  that the (number of) cells used varies between pairs. The benefit of
  this approach is that all available data is used. Use "complete.obs",
  if you want to only use the values from cells that are not `NA` in any
  of the layers. By using "masked.complete" you indicate that all layers
  have NA values in the same cells

- maxcell:

  positive integer. The maximum number of cells to be used. If this is
  smaller than ncell(x), a regular sample of `x` is used

- ...:

  additional arguments for `fun` (if it is a proper function)

## Value

If `fun` is one of the three standard statistics, you get a list with
three items: the correlation or (weighted) covariance matrix, the
(weighted) means, and the number of data cells in each comparison. The
means are also a matrix because they may depend on the combination of
layers if different cells have missing values and these are excluded
from the computation. The rows of the mean matrix represent the layer
whose (weighted) mean is being calculated and the columns represent the
layer it is being paired with. Only cells with non-missing observations
for both layers are used in the calculation of the (weighted) mean. The
diagonals of the mean and n matrices are set to missing.

If `fun` is a function, you get a single matrix.

## References

For the weighted covariance:

- Canty, M.J. and A.A. Nielsen, 2008. Automatic radiometric
  normalization of multitemporal satellite imagery with the iteratively
  re-weighted MAD transformation. Remote Sensing of Environment
  112:1025-1036.

- Nielsen, A.A., 2007. The regularized iteratively reweighted MAD method
  for change detection in multi- and hyperspectral data. IEEE
  Transactions on Image Processing 16(2):463-478.

## See also

[`global`](https://rspatial.github.io/terra/reference/global.md),
[`cov.wt`](https://rdrr.io/r/stats/cov.wt.html),
[`weighted.mean`](https://rspatial.github.io/terra/reference/weighted.mean.md)

## Examples

``` r
b <- rast(system.file("ex/logo.tif", package="terra"))   
layerCor(b, "cor")
#> $correlation
#>             red     green      blue
#> red   1.0000000 0.9980961 0.9501633
#> green 0.9980961 1.0000000 0.9658011
#> blue  0.9501633 0.9658011 1.0000000
#> 
#> $mean
#>            red    green     blue
#> red        NaN 182.2855 182.2855
#> green 185.3509      NaN 185.3509
#> blue  192.8046 192.8046      NaN
#> 
#> $n
#>        red green blue
#> red    NaN  7777 7777
#> green 7777   NaN 7777
#> blue  7777  7777  NaN
#> 

layerCor(b, "cov")
#> $covariance
#>            red    green     blue
#> red   5564.371 5443.405 4993.165
#> green 5443.405 5345.403 4974.478
#> blue  4993.165 4974.478 4962.942
#> 
#> $mean
#>            red    green     blue
#> red   182.2855 182.2855 182.2855
#> green 185.3509 185.3509 185.3509
#> blue  192.8046 192.8046 192.8046
#> 
#> $n
#>      [,1] [,2] [,3]
#> [1,] 7777 7777 7777
#> [2,] 7777 7777 7777
#> [3,] 7777 7777 7777
#> 

# weigh by column number
w <- init(b, fun="col")
layerCor(b, "weighted.cov", w=w)
#> $weighted_covariance
#>            red    green     blue
#> red   5670.750 5536.351 5009.851
#> green 5536.351 5427.161 4987.092
#> blue  5009.851 4987.092 4937.007
#> 
#> $weighted_mean
#>            red    green     blue
#> red   177.5983 177.5983 177.5983
#> green 181.3521 181.3521 181.3521
#> blue  191.5236 191.5236 191.5236
#> 

# specify another function
layerCor(b, function(x, y) cor(x, y, method="spearman"))
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.9884613 0.9291036
#> [2,] 0.9884613 1.0000000 0.9425167
#> [3,] 0.9291036 0.9425167 1.0000000
```

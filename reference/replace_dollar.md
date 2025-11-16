# Replace with \$\<-

Replace a layer of a SpatRaster, or an attribute variable of a
SpatVector

## Usage

``` r
# S4 method for class 'SpatRaster'
x$name <- value

# S4 method for class 'SpatVector'
x$name <- value

# S4 method for class 'SpatExtent'
x$name <- value
```

## Arguments

- x:

  SpatRaster, SpatVector or SpatExtent

- name:

  character. If `x` is a SpatRaster: layer name. If `x` is a SpatVector:
  variable name. If `x` is a SpatExtent: "xmin", "xmax". "ymin" or
  "ymax"

- value:

  if `x` is a SpatRaster, a SpatRaster for which this `TRUE`:
  `nlyr(value) == length(i)`; if `x` is a SpatVector, a vector of new
  values; if `x` is a SpatExtent a single number

## Value

Same as `x`

## See also

` [[<-, [<-, `[`$`](https://rspatial.github.io/terra/reference/subset_dollar.md)

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
v$ID_1 <- LETTERS[1:12]
v$new <- sample(12)
values(v)
#>    ID_1       NAME_1 ID_2           NAME_2 AREA    POP new
#> 1     A     Diekirch    1         Clervaux  312  18081  10
#> 2     B     Diekirch    2         Diekirch  218  32543  11
#> 3     C     Diekirch    3          Redange  259  18664   7
#> 4     D     Diekirch    4          Vianden   76   5163   5
#> 5     E     Diekirch    5            Wiltz  263  16735  12
#> 6     F Grevenmacher    6       Echternach  188  18899   1
#> 7     G Grevenmacher    7           Remich  129  22366   4
#> 8     H Grevenmacher   12     Grevenmacher  210  29828   8
#> 9     I   Luxembourg    8         Capellen  185  48187   6
#> 10    J   Luxembourg    9 Esch-sur-Alzette  251 176820   2
#> 11    K   Luxembourg   10       Luxembourg  237 182607   3
#> 12    L   Luxembourg   11           Mersch  233  32112   9
```

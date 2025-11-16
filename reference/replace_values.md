# Replace values of a SpatRaster

Replace values of a SpatRaster. These are convenience functions for
smaller objects only. For larger rasters see `link{classify}` or
[`subst`](https://rspatial.github.io/terra/reference/subst.md)

## Usage

``` r
# S4 method for class 'SpatRaster,ANY,ANY,ANY'
x[i, j, k] <- value

# S4 method for class 'SpatVector,ANY,ANY'
x[i, j] <- value

# S4 method for class 'SpatExtent,numeric,missing'
x[i, j] <- value
```

## Arguments

- x:

  SpatRaster

- i:

  row numbers. numeric, logical, or missing for all rows. Can also be a
  SpatRaster or SpatVector

- j:

  column numbers. numeric, logical or missing for all columns

- k:

  layer number. numeric, logical or missing for all layers

- value:

  numeric, matrix, or data.frame

## Value

SpatRaster

## See also

[`classify`](https://rspatial.github.io/terra/reference/classify.md)`, `[`subst`](https://rspatial.github.io/terra/reference/subst.md)`, `[`set.values`](https://rspatial.github.io/terra/reference/inplace.md)`, `[`values`](https://rspatial.github.io/terra/reference/values.md)`, [[<-`

## Examples

``` r
## SpatRaster
r <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
r[] <- 1:25
r[1,] <- 5
r[,2] <- 10
r[r>10] <- NA

## SpatVector
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
v[2,2] <- "hello"
v[1,] <- v[10,]
v[,3] <- v[,1]
v[2, "NAME_2"] <- "terra"
head(v, 3)
#>   ID_1     NAME_1 ID_2           NAME_2 AREA    POP
#> 1    3 Luxembourg    3 Esch-sur-Alzette  251 176820
#> 2    1      hello    1            terra  218  32543
#> 3    1   Diekirch    1          Redange  259  18664
```

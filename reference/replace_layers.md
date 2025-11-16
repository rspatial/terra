# Replace layers or variables

Replace the layers of SpatRaster with (layers from) another SpatRaster
or replace variables of a SpatVector. You can also create new
layers/variables with these methods.

## Usage

``` r
# S4 method for class 'SpatRaster,numeric'
x[[i]] <- value

# S4 method for class 'SpatRaster,character'
x[[i]] <- value

# S4 method for class 'SpatVector,numeric'
x[[i]] <- value

# S4 method for class 'SpatVector,character'
x[[i]] <- value
```

## Value

SpatRaster

## Arguments

- x:

  SpatRaster or SpatVector

- i:

  if `x` is a SpatRaster: layer number(s) of name(s). If `x` is a
  SpatVector: variable number(s) or name(s) (column of the attributes)

- value:

  if `x` is a SpatRaster: SpatRaster for which this `TRUE`:
  `nlyr(value) == length(i)`. if `x` is a SpatVector: vector or
  data.frame

## See also

`$<-, [<-`

## Examples

``` r
# raster
s <- rast(system.file("ex/logo.tif", package="terra"))   
s[["red"]] <- mean(s)
s[[2]] <- sqrt(s[[1]])

# vector
v <- vect(system.file("ex/lux.shp", package="terra"))
v[["ID_1"]] <- 12:1
```

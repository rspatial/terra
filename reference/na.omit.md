# Find and remove geometries that are NA

Find geometries that are NA; or remove geometries and/or records that
are `NA`.

## Usage

``` r
# S4 method for class 'SpatVector'
is.na(x)

# S4 method for class 'SpatVector'
na.omit(object, field=NA, geom=FALSE)
```

## Arguments

- x:

  SpatVector

- object:

  SpatVector

- field:

  character or NA. If `NA`, missing values in the attributes are
  ignored. Other values are either one or more field (variable) names,
  or `""` to consider all fields

- geom:

  logical. If `TRUE` empty geometries are removed

## Value

SpatVector

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
v$test <- c(1,2,NA)
nrow(v)
#> [1] 12
x <- na.omit(v, "test")
nrow(x)
#> [1] 8
```

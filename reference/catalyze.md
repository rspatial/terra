# Factors to numeric

Change a categorical layer into one or more numerical layers. With
`as.numeric` you can transfer the active category values to cell values
in a non-categorical SpatRaster. `catalyze` creates new layers for each
category.

## Usage

``` r
# S4 method for class 'SpatRaster'
as.numeric(x, index=NULL, filename="", ...)

# S4 method for class 'SpatRaster'
catalyze(x, filename="", ...)
```

## Arguments

- x:

  SpatRaster

- index:

  positive integer or category indicating the category to use. If `NULL`
  the active category is used

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## Author

Andrew Gene Brown, Robert J. Hijmans

## See also

[`activeCat`](https://rspatial.github.io/terra/reference/activeCat.md),
[`cats`](https://rspatial.github.io/terra/reference/factors.md)

## Examples

``` r
set.seed(0)
r <- rast(nrows=10, ncols=10)
values(r) <- sample(3, ncell(r), replace=TRUE) + 10
d <- data.frame(id=11:13, cover=c("forest", "water", "urban"), letters=letters[1:3], value=10:12)
levels(r) <- d
catalyze(r)
#> class       : SpatRaster 
#> size        : 10, 10, 3  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> names       : cover, letters, value 
#> min values  :     1,       1,    10 
#> max values  :     3,       3,    12 

activeCat(r) <- 3
as.numeric(r)
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> name        : letters 
#> min value   :       1 
#> max value   :       3 
```

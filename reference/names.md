# Names of Spat\* objects

Get or set the names of the layers of a SpatRaster or the attributes of
a SpatVector.

See [`set.names`](https://rspatial.github.io/terra/reference/inplace.md)
for in-place setting of names.

## Usage

``` r
# S4 method for class 'SpatRaster'
names(x)

# S4 method for class 'SpatRaster'
names(x) <- value

# S4 method for class 'SpatRasterDataset'
names(x)

# S4 method for class 'SpatRasterDataset'
names(x) <- value

# S4 method for class 'SpatVector'
names(x)

# S4 method for class 'SpatVector'
names(x) <- value
```

## Arguments

- x:

  SpatRaster, SpatRasterDataset, or SpatVector

- value:

  character (vector)

## Value

character

## Note

terra enforces neither unique nor valid names. See
[`make.unique`](https://rdrr.io/r/base/make.unique.html) to create
unique names and [`make.names`](https://rdrr.io/r/base/make.names.html)
to make syntactically valid names.

## Examples

``` r
s <- rast(ncols=5, nrows=5, nlyrs=3)
nlyr(s)
#> [1] 3
names(s)
#> [1] "lyr.1" "lyr.2" "lyr.3"
names(s) <- c("a", "b", "c")
names(s)
#> [1] "a" "b" "c"

# SpatVector names
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
names(v)
#> [1] "ID_1"   "NAME_1" "ID_2"   "NAME_2" "AREA"   "POP"   
names(v) <- paste0(substr(names(v), 1, 2), "_", 1:ncol(v))
names(v)
#> [1] "ID_1" "NA_2" "ID_3" "NA_4" "AR_5" "PO_6"
```

# Read all cell values into memory

Reads all cell values of a SpatRaster or SpatRasterDataset into memory.

Using this method is discouraged as it is not necessary for processing
the data and may lead to excessive memory use that will slow down your
computer or worse. It cannot be used for SpatRasters that are based on
very large files.

The method may be useful if a relatively small dataset is used
repeatedly, such that efficiency gains are made because the values only
need to be read from disk once.

## Usage

``` r
# S4 method for class 'SpatRaster'
toMemory(x)

# S4 method for class 'SpatRasterDataset'
toMemory(x)
```

## Arguments

- x:

  SpatRaster or SpatRasterDataset

## Value

Same as `x`

## See also

[`values`](https://rspatial.github.io/terra/reference/values.md)`, `[`as.data.frame`](https://rspatial.github.io/terra/reference/as.data.frame.md)`, `[`readValues`](https://rspatial.github.io/terra/reference/readwrite.md)`, `[`inMemory`](https://rspatial.github.io/terra/reference/sources.md)

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
sources(r)
#> [1] "/home/runner/work/_temp/Library/terra/ex/elev.tif"
inMemory(r)
#> [1] FALSE
x <- toMemory(r)
inMemory(x)
#> [1] TRUE
```

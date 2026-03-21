# Set the NA flag

The main purpose of this method is to allow correct reading of a
SpatRaster that is based on a file that has an incorrect or missing NA
flag. The file is not changed, but the flagged value is set to NA when
values are read from the file. In contrast, if the cell values are in
memory the change is made immediately.

To change cell values, you can use
[`classify`](https://rspatial.github.io/terra/reference/classify.md)

## Usage

``` r
# S4 method for class 'SpatRaster'
NAflag(x)

# S4 method for class 'SpatRaster'
NAflag(x) <- value
```

## Arguments

- x:

  SpatRaster

- value:

  numeric. The value to be interpreted as NA; set this before reading
  the values from the file. This can be a single value, or multiple
  values, one for each data source (file / subdataset)

## Value

none or numeric

## See also

[`classify`](https://rspatial.github.io/terra/reference/classify.md)

## Examples

``` r
s <- rast(system.file("ex/logo.tif", package="terra"))[[1]]   
NAflag(s) <- 255
plot(s)

NAflag(s)
#> [1] 255
```

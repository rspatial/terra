# segregate

Create a SpatRaster with a layer for each class (value, or subset of the
values) in the input SpatRaster. For example, if the input has
vegetation types, this function will create a layer (presence/absence;
dummy variable) for each of these classes.

This is called "one-hot encoding" or "dummy encoding" (for a dummy
encoding scheme you can remove (any) one of the output layers as it is
redundant).

## Usage

``` r
# S4 method for class 'SpatRaster'
segregate(x, classes=NULL, keep=FALSE, other=0, round=FALSE, digits=0, filename="", ...)
```

## Arguments

- x:

  SpatRaster

- classes:

  numeric. The values (classes) for which layers should be made. If
  `NULL` all classes are used

- keep:

  logical. If `TRUE`, cells that are of the class represented by a layer
  get that value, rather than a value of 1

- other:

  numeric. Value to assign to cells that are not of the class
  represented by a layer

- round:

  logical. Should the values be rounded first?

- digits:

  integer. Number of digits to round the values to

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`split`](https://rspatial.github.io/terra/reference/split.md)

## Examples

``` r
r <- rast(nrows=5, ncols=5)
values(r) <- rep(c(1:4, NA), each=5)
b <- segregate(r)
bb <- segregate(r, keep=TRUE, other=NA)
```

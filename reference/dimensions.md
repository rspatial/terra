# Dimensions of a SpatRaster or SpatVector and related objects

Get the number of rows (`nrow`), columns (`ncol`), cells (`ncell`),
layers (`nlyr`), sources (`nsrc`), the size `size` (`nlyr(x)*ncell(x)`),
or spatial resolution of a SpatRaster.

`length` returns the number of sub-datasets in a SpatRasterDataset or
SpatVectorCollection.

For a SpatVector `length(x)` is the same as `nrow(x)`.

You can also set the number of rows or columns or layers. When setting
dimensions, all cell values are dropped.

## Usage

``` r
# S4 method for class 'SpatRaster'
ncol(x)

# S4 method for class 'SpatRaster'
nrow(x)

# S4 method for class 'SpatRaster'
nlyr(x)

# S4 method for class 'SpatRaster'
ncell(x)

# S4 method for class 'SpatRaster'
nsrc(x)

# S4 method for class 'SpatRaster,numeric'
ncol(x) <- value

# S4 method for class 'SpatRaster,numeric'
nrow(x) <- value

# S4 method for class 'SpatRaster,numeric'
nlyr(x) <- value

# S4 method for class 'SpatRaster'
res(x)

# S4 method for class 'SpatRaster,numeric'
res(x) <- value

# S4 method for class 'SpatRaster'
xres(x)

# S4 method for class 'SpatRaster'
yres(x)

# S4 method for class 'SpatVector'
ncol(x)

# S4 method for class 'SpatVector'
nrow(x)

# S4 method for class 'SpatVector'
length(x)
```

## Arguments

- x:

  SpatRaster or SpatVector or related objects

- value:

  For ncol and nrow: positive integer. For res: one or two positive
  numbers

## Value

integer

## See also

[ext](https://rspatial.github.io/terra/reference/ext.md)

## Examples

``` r
r <- rast()
ncol(r)
#> [1] 360
nrow(r)
#> [1] 180
nlyr(r)
#> [1] 1
dim(r)
#> [1] 180 360   1
nsrc(r)
#> [1] 1
ncell(r)
#> [1] 64800

rr  <- c(r,r)
nlyr(rr)
#> [1] 2
nsrc(rr)
#> [1] 2
ncell(rr)
#> [1] 64800

nrow(r) <- 18
ncol(r) <- 36
# equivalent to
dim(r) <- c(18, 36) 

dim(r)
#> [1] 18 36  1
dim(r) <- c(10, 10, 5)
dim(r)
#> [1] 10 10  5

xres(r)
#> [1] 36
yres(r)
#> [1] 18
res(r)
#> [1] 36 18

res(r) <- 1/120
# different xres and yres
res(r) <- c(1/120, 1/60)
```

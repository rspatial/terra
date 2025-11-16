# select cells with high or low values

Identify n cells that have the highest or lowest values in the first
layer of a SpatRaster.

## Usage

``` r
# S4 method for class 'SpatRaster'
selectHighest(x, n, low=FALSE)
```

## Arguments

- x:

  SpatRaster. Only the first layer is processed

- n:

  The number of cells to select

- low:

  logical. If `TRUE`, the lowest values are selected instead of the
  highest values

## Value

SpatRaster

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
x <- selectHighest(r, 1000)
y <- selectHighest(r, 1000, TRUE)

m <- merge(y-1, x)
levels(m) <- data.frame(id=0:1, elevation=c("low", "high"))
plot(m)
```

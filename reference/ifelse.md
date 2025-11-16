# ifelse for SpatRasters

Implementation of [`ifelse`](https://rdrr.io/r/base/ifelse.html) for
SpatRasters. This method allows for a concise expression of what can
otherwise be achieved with a combination of
[`classify`](https://rspatial.github.io/terra/reference/classify.md),
[`mask`](https://rspatial.github.io/terra/reference/mask.md), and
[`cover`](https://rspatial.github.io/terra/reference/cover.md).

`ifel` is an `R` equivalent to the `Con` method in ArcGIS (arcpy).

## Usage

``` r
# S4 method for class 'SpatRaster'
ifel(test, yes, no, filename="", ...)
```

## Arguments

- test:

  SpatRaster with logical (TRUE/FALSE) values

- yes:

  SpatRaster or numeric

- no:

  SpatRaster or numeric

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## Examples

``` r
r <- rast(nrows=5, ncols=5, xmin=0, xmax=1, ymin=0, ymax=1)
values(r) <- c(-10:0, NA, NA, NA, 0:10)

x <- ifel(r > 1, 1, r)
# same as 
a <- classify(r, cbind(1, Inf, 1))
# or
b <- app(r, fun=function(i) {i[i > 1] <- 1; i})
# or 
d <- clamp(r, -Inf, 1)
# or (not recommended for large datasets)
e <- r
e[e>1] <- 1

## other examples
f <- ifel(is.na(r), 100, r)

z <- ifel(r > -2 & r < 2, 100, 0)

# nested expressions
y <- ifel(r > 1, 1, ifel(r < -1, -1, r))

k <- ifel(r > 0, r+10, ifel(r < 0, r-10, 3))
```

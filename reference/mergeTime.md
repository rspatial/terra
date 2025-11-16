# merge SpatRasters by timelines to create a single timeseries

Combine SpatRasters with partly overlapping time-stamps to create a
single time series. If there is no overlap between the SpatRasters there
is no point in using this function (use
[`c`](https://rspatial.github.io/terra/reference/c.md) instead).

Also note that time gaps are not filled. You can use
[`fillTime`](https://rspatial.github.io/terra/reference/fillTime.md) to
do that.

## Usage

``` r
# S4 method for class 'SpatRasterDataset'
mergeTime(x, fun=mean, filename="", ...)
```

## Arguments

- x:

  SpatRasterDataset

- fun:

  A function that reduces a vector to a single number, such as `mean` or
  `min`

- filename:

  character. Output filename

- ...:

  list with named options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## Examples

``` r
r <- rast(system.file("ex/logo.tif", package="terra"))   
s1 <- c(r, r)
time(s1) <- as.Date("2001-01-01") + 0:5
s1 <- s1/10
time(s1) <- as.Date("2001-01-07") + 0:5
s2 <- s1*10
time(s2) <- as.Date("2001-01-05") + 0:5
x <- sds(s1, s1, s2)

m <- mergeTime(x, mean)
```

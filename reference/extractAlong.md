# extract values along lines

Extract raster values along a line. That is, the returned values are
ordered along the line. That is not the case with
[`extract`](https://rspatial.github.io/terra/reference/extract.md)

## Usage

``` r
extractAlong(x, y, ID=TRUE, cells=FALSE, xy=FALSE, online=FALSE, bilinear=TRUE)
```

## Arguments

- x:

  SpatRaster

- y:

  SpatVector with lines geometry

- ID:

  logical. Should an ID column be added? If so, the first column
  returned has the IDs (record numbers) of input SpatVector `y`

- cells:

  logical. If `TRUE` the cell numbers are also returned

- xy:

  logical. If `TRUE` the coordinates of the cells traversed by `y` are
  also returned. See
  [`xyFromCell`](https://rspatial.github.io/terra/reference/xyCellFrom.md)

- online:

  logical. If `TRUE` the returned coordinates are snapped to `y`

- bilinear:

  logical. If `TRUE` the returned raster values computed with bilinear
  interpolation from the nearest four cells. Only relevant if
  `online=TRUE`

## Value

data.frame

## See also

[`extract`](https://rspatial.github.io/terra/reference/extract.md)

## Examples

``` r
r <- rast(ncols=36, nrows=18, vals=1:(18*36))
cds1 <- rbind(c(-50,0), c(0,60), c(40,5), c(15,-45), c(-10,-25))
cds2 <- rbind(c(80,20), c(140,60), c(160,0), c(140,-55))
lines <- vect(list(cds1, cds2), "lines")

extractAlong(r, lines)
#>    ID lyr.1
#> 1   1   337
#> 2   1   302
#> 3   1   266
#> 4   1   267
#> 5   1   231
#> 6   1   232
#> 7   1   196
#> 8   1   197
#> 9   1   161
#> 10  1   162
#> 11  1   126
#> 12  1    91
#> 13  1   127
#> 14  1   163
#> 15  1   164
#> 16  1   200
#> 17  1   201
#> 18  1   237
#> 19  1   273
#> 20  1   274
#> 21  1   310
#> 22  1   310
#> 23  1   346
#> 24  1   382
#> 25  1   381
#> 26  1   417
#> 27  1   453
#> 28  1   452
#> 29  1   488
#> 30  1   488
#> 31  1   487
#> 32  1   451
#> 33  1   450
#> 34  1   414
#> 35  2   279
#> 36  2   243
#> 37  2   244
#> 38  2   208
#> 39  2   209
#> 40  2   210
#> 41  2   174
#> 42  2   175
#> 43  2   139
#> 44  2   140
#> 45  2   105
#> 46  2   141
#> 47  2   177
#> 48  2   213
#> 49  2   214
#> 50  2   250
#> 51  2   286
#> 52  2   322
#> 53  2   323
#> 54  2   358
#> 55  2   394
#> 56  2   430
#> 57  2   429
#> 58  2   465
#> 59  2   501
#> 60  2   537
```

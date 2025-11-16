# k_means

Compute k-means clusters for a SpatRaster. For large SpatRasters (with
`ncell(x) > maxcell`) this is done in two steps. First a sample of the
cells is used to compute the cluster centers. Then each cell is assigned
to a cluster by computing the distance to these centers.

## Usage

``` r
# S4 method for class 'SpatRaster'
k_means(x, centers=3, ..., maxcell=1000000, filename="", overwrite=FALSE, wopt=list())
```

## Arguments

- x:

  SpatRaster

- centers:

  either the number of clusters, or a set of initial (distinct) cluster
  centres. If a number, a random set of (distinct) cells in `x` is
  chosen as the initial centres

- ...:

  additional arguments passed to
  [`kmeans`](https://rdrr.io/r/stats/kmeans.html)

- maxcell:

  positive integer. The size of the regular sample used if it is smaller
  than `ncell(x)`

- filename:

  character. Output filename (ignored if `as.raster=FALSE`)

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- wopt:

  list with additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`kmeans`](https://rdrr.io/r/stats/kmeans.html)

## Examples

``` r
f <- system.file("ex/logo.tif", package = "terra")
r <- rast(f)
km <- k_means(r, centers=5)
km
#> class       : SpatRaster 
#> size        : 77, 101, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source(s)   : memory
#> name        : lyr1 
#> min value   :    1 
#> max value   :    5 
```

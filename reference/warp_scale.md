# Compute warp resampling scale

Compute the ratio between source and destination pixel sizes for use
with [`project`](https://rspatial.github.io/terra/reference/project.md).
These values correspond to the GDAL warp options XSCALE and YSCALE. By
default, GDAL computes them independently for each processing chunk,
which can produce visible block-boundary artifacts when the ratio varies
across the raster (particularly in projections with significant
deformation). Setting fixed values via `project(..., xscale=, yscale=)`
eliminates these artifacts.

The function samples a grid of points in the source raster, projects
them to the destination CRS, and returns summary statistics of the local
resampling ratios.

## Usage

``` r
warp_scale(x, y, n=21)
```

## Arguments

- x:

  SpatRaster. The source raster

- y:

  SpatRaster or character. A template raster or a CRS description for
  the destination

- n:

  integer. Number of sample points along each axis (the total number of
  sample points is `n*n`)

## Value

A list with two named numeric vectors (`xscale` and `yscale`), each
containing the quantiles (0%, 25%, 50%, 75%, 100%) of the local
resampling ratio across the sampled points. The median (`"50%"`) is
typically the best choice for a global fixed scale. Values equal one for
no resampling, below one for downsampling, and above one for upsampling.

## See also

[`project`](https://rspatial.github.io/terra/reference/project.md)

## Examples

``` r
if (FALSE) { # \dontrun{
a <- rast(ncols=360, nrows=180, xmin=-180, xmax=180, ymin=-90, ymax=90,
          crs="+proj=longlat +datum=WGS84")
values(a) <- 1:ncell(a)

sc <- warp_scale(a, "+proj=robin")
sc
# Use the median values
b <- project(a, "+proj=robin", xscale=sc$xscale["50%"], yscale=sc$yscale["50%"])
} # }
```

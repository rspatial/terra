# Distance on a grid

The function calculates the distance to cells of a SpatRaster when the
path has to go through the centers of the eight neighboring raster
cells.

The default distance (when `scale=1`, is meters if the coordinate
reference system (CRS) of the SpatRaster is longitude/latitude
(`+proj=longlat`) and in the linear units of the CRS (typically meters)
in other cases.

Distances are computed by summing local distances between cells, which
are connected with their neighbors in 8 directions.

The shortest distance to the cells with the `target` value is computed
for all cells that are not `NA`. Cells that are `NA` cannot be traversed
and are ignored, unless the target itself is `NA`, in which case the
distance to the nearest cell that is not `NA` is computed for all cells
that are `NA`.

## Usage

``` r
# S4 method for class 'SpatRaster'
gridDist(x, target=0, scale=1, maxiter=50, filename="", ...)
```

## Arguments

- x:

  SpatRaster

- target:

  numeric. value of the target cells (where to compute distance to)

- scale:

  numeric. Scale factor. For longitude/latitude data 1 = "m" and 1000 =
  "km". For planar data that is also the case of the distance unit of
  the crs is "m"

- maxiter:

  numeric. The maximum number of iterations. Increase this number if you
  get the warning that `costDistance` did not converge. Only relevant
  when target is not `NA`

- filename:

  character. output filename (optional)

- ...:

  additional arguments as for
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## See also

See [`distance`](https://rspatial.github.io/terra/reference/distance.md)
for "as the crow flies" distance, and
[`costDist`](https://rspatial.github.io/terra/reference/costDist.md) for
distance across a landscape with variable friction

## Value

SpatRaster

## Examples

``` r
# global lon/lat raster
r <- rast(ncol=10,nrow=10, vals=1)
r[48] <- 0
r[66:68] <- NA
d <- gridDist(r) 
plot(d)



# planar
crs(r) <- "+proj=utm +zone=15 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
d <- gridDist(r) 
plot(d)


# distance to cells that are not NA 
rr <- classify(r, cbind(1, NA))
dd <- gridDist(rr, NA) 

```

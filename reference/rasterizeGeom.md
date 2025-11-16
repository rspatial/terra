# Rasterize geometric properties of vector data

Rasterization of geometric properties of vector data. You can get the
count of the number of geometries in each cell; the area covered by
polygons; the length of the lines; or the number of lines that cross the
boundary of each cell. See
[`rasterize`](https://rspatial.github.io/terra/reference/rasterize.md)
for standard rasterization (of attribute values associated with
geometries).

The area of polygons is intended for summing the area of polygons that
are relatively small relative to the raster cells, and for when there
may be multiple polygons per cell. See `rasterize(fun="sum")` for
counting large polygons and `rasterize(cover=TRUE)` to get the fraction
that is covered by larger polygons.

## Usage

``` r
# S4 method for class 'SpatVector,SpatRaster'
rasterizeGeom(x, y, fun="count", unit="m", filename="", ...)
```

## Arguments

- x:

  SpatVector

- y:

  SpatRaster

- fun:

  character. "count", "area", "length", or "crosses"

- unit:

  character. "m" or "km"

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## See also

[`rasterize`](https://rspatial.github.io/terra/reference/rasterize.md)

## Value

SpatRaster

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
r <- rast(v, res=.1)

# length of lines
lns <- as.lines(v)
x <- rasterizeGeom(lns, r, fun="length", "km")

# count of points
set.seed(44)
pts <- spatSample(v, 100)
y <- rasterizeGeom(pts, r)

# area of polygons
pols <- buffer(pts, 1000)
z <- rasterizeGeom(pols, r, fun="area")
```

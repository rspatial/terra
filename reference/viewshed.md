# Compute a viewshed

Use elevation data to compute the locations that can be seen, or how
much higher they would have to be to be seen, from a certain position.
The raster data coordinate reference system must be planar (not
lon/lat), with the elevation values in the same unit as the distance
unit of the coordinate reference system.

## Usage

``` r
# S4 method for class 'SpatRaster'
viewshed(x, loc, observer=1.80, target=0, curvcoef=6/7, output="yes/no", filename="", ...)
```

## Arguments

- x:

  SpatRaster, single layer with elevation values. Values should have the
  same unit as the map units

- loc:

  location (x and y coordinates) or a cell number

- observer:

  numeric. The height above the elevation data of the observer

- target:

  numeric. The height above the elevation data of the targets

- curvcoef:

  numeric. Coefficient to consider the effect of the curvature of the
  earth and refraction of the atmosphere. The elevation values are
  corrected with:
  `elevation = elevation - curvcoeff * (distance)^2 / (earth_diameter)`.
  This means that with the default value of 0.85714, you lose sight of
  about 1 meter of elevation for each 385 m of planar distance

- output:

  character. Can be "yes/no" to get a binary (logical) output showing
  what areas are visible; "land" to get the height above the current
  elevation that would be visible; or "sea" the elevation above sea
  level that would be visible

- filename:

  character. Output filename

- ...:

  Options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## See also

[`terrain`](https://rspatial.github.io/terra/reference/terrain.md)

## References

The algorithm used is by Wang et al.:
https://www.asprs.org/wp-content/uploads/pers/2000journal/january/2000_jan_87-90.pdf.

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
x <- project(r, "EPSG:2169")
p <- cbind(70300, 96982)
v <- viewshed(x, p, 0, 0, 0.85714)
```

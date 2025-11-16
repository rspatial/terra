# Compute surface area from elevation data

It is often said that if Wales was flattened out it would have an area
bigger than England. This function computes the surface area for a
raster with elevation values, taking into account the sloping nature of
the surface.

## Usage

``` r
# S4 method for class 'SpatRaster'
surfArea(x, filename="", ...)
```

## Arguments

- x:

  SpatRaster with elevation values. Currently the raster CRS must be
  planar and have the same distance units (e.g. m) as the elevation
  values

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`expanse`](https://rspatial.github.io/terra/reference/expanse.md),
[`cellSize`](https://rspatial.github.io/terra/reference/cellSize.md)

## References

Jenness, Jeff S., 2004. Calculating Landscape Surface Area from Digital
Elevation Models. Wildlife Society Bulletin 32(3): 829-839

## Author

Barry Rowlingson

## Examples

``` r
v <- rast(volcano, crs="local")
x <- terra::surfArea(v)
```

# Rotate data along longitude

Rotate a SpatRaster that has longitude coordinates from 0 to 360, to
standard coordinates between -180 and 180 degrees (or vice-versa).
Longitude between 0 and 360 is frequently used in global climate models.

Rotate a SpatVector as for a SpatRaster with, or with `split=FALSE` to
correct for coordinates that are connected across the date line (and end
up at the "other side" of the longitude scale).

## Usage

``` r
# S4 method for class 'SpatRaster'
rotate(x, filename="", ...)

# S4 method for class 'SpatVector'
rotate(x, longitude=0, split=TRUE, left=TRUE, normalize=FALSE)
```

## Arguments

- x:

  SpatRaster or SpatVector

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

- longitude:

  numeric. The longitude around which to rotate

- split:

  logical. Should geometries be split at `longitude`?

- left:

  logical. Rotate to the left or to the right?

- normalize:

  logical. Should the output be normalized to longitudes between -180
  and 180? See
  [`normalize.longitude`](https://rspatial.github.io/terra/reference/normalize.longitude.md)

## Value

SpatRaster

## See also

[`shift`](https://rspatial.github.io/terra/reference/shift.md) and
[`spin`](https://rspatial.github.io/terra/reference/spin.md)

## Examples

``` r
x <- rast(nrows=9, ncols=18, nl=3, xmin=0, xmax=360)
v <- rep(as.vector(t(matrix(1:ncell(x), nrow=9, ncol=18))), 3)
values(x) <- v
z <- rotate(x)

if (FALSE) { # \dontrun{
#SpatVector
p <- rbind(c(3847903, 1983584 ), c(3847903, 5801864), c(8301883, 5801864), c(8301883, 1983584 ))
p <- vect(p, "polygons", crs="+init=EPSG:3347")
d <- densify(p, 100000)
g <- project(d, "+proj=longlat")

x <- rotate(g, 50)
plot(g)
lines(x, col="red")
} # }

## rotate countries to 0-360 longitude
#w <- geodata::world(path=".")
#x <- rotate(w, long=0, split=TRUE, left=FALSE)
```

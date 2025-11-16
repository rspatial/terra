# Zoom in on a map

Zoom in on a map (plot) by providing a new extent, by default this is
done by clicking twice on the map.

## Usage

``` r
# S4 method for class 'SpatRaster'
zoom(x, e=draw(), maxcell=100000, layer=1, new=FALSE, ...)

# S4 method for class 'SpatVector'
zoom(x, e=draw(), new=FALSE, ...)
```

## Arguments

- x:

  SpatRaster

- e:

  SpatExtent

- maxcell:

  positive integer. Maximum number of cells used for the map

- layer:

  positive integer to select the layer to be used

- new:

  logical. If `TRUE`, the zoomed in map will appear on a new device
  (window)

- ...:

  additional arguments passed to
  [`plot`](https://rspatial.github.io/terra/reference/plot.md)

## Value

SpatExtent (invisibly)

## See also

[`draw`](https://rspatial.github.io/terra/reference/draw.md),
[`plot`](https://rspatial.github.io/terra/reference/plot.md)

# Catchment delineation

delineate the area covered by a catchment from a SpatRaster with flow
direction and a pour-point (catchment outlet).

## Usage

``` r
# S4 method for class 'SpatRaster'
watershed(x, pourpoint, filename="",...)
```

## Arguments

- x:

  SpatRaster with flow direction. See
  [`terrain`](https://rspatial.github.io/terra/reference/terrain.md).

- pourpoint:

  matrix or SpatVector with the pour point location

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## Author

Ezio Crestaz, Emanuele Cordano, Roman Seliger

## Examples

``` r
elev <- rast(system.file('ex/elev_vinschgau.tif', package="terra"))
flowdir <- terrain(elev, "flowdir")
## pour point at Naturns 
pp <- cbind(653358.3, 5168222)
w <- watershed(flowdir, pp)
```

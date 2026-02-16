# Create a graticule

Create a graticule. That is, a grid of lon/lat lines that can be used on
a projected map.

The object returned, a SpatGraticule, can be plotted with `plot` and
`lines`. There is also a `crop` method.

## Usage

``` r
graticule(lon=30, lat=30, crs="")
```

## Arguments

- lon:

  numeric. Either a single number (the interval between longitudes), or
  a vector with longitudes

- lat:

  numeric. Either a single number (the interval between latitudes), or a
  vector with latitudes

- crs:

  character. The coordinate reference system to use

## Value

SpatGraticule

## See also

[`plot<SpatGraticule>`](https://rspatial.github.io/terra/reference/plot_graticule.md).

## Examples

``` r
g <- graticule(60, 30, crs="+proj=robin")
g
#> class       : SpatGraticule 
#> lon         : -180 -120 -60 0 60 120 180 
#> lat         : -90 -60 -30 0 30 60 90 
#> coord. ref. : +proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs 
#> extent      : -17005833, 17005833, -8625155, 8625155  (xmin, xmax, ymin, ymax)

graticule(90, c(-90, -60, -23.5, 0, 23.5, 60, 90), crs="+proj=robin")
#> class       : SpatGraticule 
#> lon         : -180 -90 0 90 180 
#> lat         : -90 -60 -23.5 0 23.5 60 90 
#> coord. ref. : +proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs 
#> extent      : -17005833, 17005833, -8625155, 8625155  (xmin, xmax, ymin, ymax)
```

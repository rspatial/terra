# Class "SpatExtent"

Objects of class SpatExtent are used to define the spatial extent
(extremes) of objects of the SpatRaster class.

## Objects from the Class

You can use the
[`ext`](https://rspatial.github.io/terra/reference/ext.md) function to
create SpatExtent objects, or to extract them from a SpatRaster,
SpatVector or related objects.

## Methods

- show:

  display values of a SpatExtent object

## Examples

``` r
e <- ext(-180, 180, -90, 90)
e
#> SpatExtent : -180, 180, -90, 90 (xmin, xmax, ymin, ymax)
```

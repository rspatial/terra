# Extract values for a range of layers from a SpatRaster

Extract values from a SpatRaster for a set of locations and a range of
layers. To extract values for a single or all layers, use
[`extract`](https://rspatial.github.io/terra/reference/extract.md)

## Usage

``` r
# S4 method for class 'SpatRaster'
extractRange(x, y, first, last, lyr_fun=NULL, 
    geom_fun=NULL, ID=FALSE, na.rm=TRUE, bind=FALSE, ...)
```

## Arguments

- x:

  SpatRaster

- y:

  SpatVector (points, lines, or polygons). Alternatively, for points, a
  2-column matrix or data.frame (x, y) or (lon, lat). Or a vector with
  cell numbers

- first:

  layer name of number, indicating the first layer in the range of
  layers to be considered

- last:

  layer name or number, indicating the last layer in the range to be
  considered

- lyr_fun:

  function to summarize the extracted data across layers

- geom_fun:

  function to summarize the extracted data for each line or polygon
  geometry. Ignored if `y` has point geometry

- ID:

  logical. Should an ID column be added? If so, the first column
  returned has the IDs (record numbers) of `y`

- na.rm:

  logical. Should missing values be ignored?

- bind:

  logical. If `TRUE`, the extracted values are `cbind`-ed to `y`

- ...:

  additional arguments passed to `extract`

## Value

numeric or data.frame

## See also

[`extract`](https://rspatial.github.io/terra/reference/extract.md)

## Examples

``` r
r <- rast(system.file("ex/logo.tif", package="terra"))   
xy <- data.frame(lon=c(50,80), lat=c(30, 60))
extract(r, xy)
#> Warning: Cannot find coordinate operations from `GEOGCRS["unknown",DATUM["World Geodetic System 1984",ELLIPSOID["WGS 84",6378137,298.257223563,LENGTHUNIT["metre",1]],ID["EPSG",6326]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]],CS[ellipsoidal,2],AXIS["longitude",east,ORDER[1],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]],AXIS["latitude",north,ORDER[2],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]]' to `ENGCRS["Cartesian (Meter)",EDATUM["Unknown engineering datum"],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]' (GDAL error 6)
#> Warning: [extract] transforming vector data to the CRS of the raster
#> Warning: number of rows of result is not a multiple of vector length (arg 1)
#> [1] ID    red   green blue 
#> <0 rows> (or 0-length row.names)
extract(r, xy, layer=c("red", "green"))
#> Warning: Cannot find coordinate operations from `GEOGCRS["unknown",DATUM["World Geodetic System 1984",ELLIPSOID["WGS 84",6378137,298.257223563,LENGTHUNIT["metre",1]],ID["EPSG",6326]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]],CS[ellipsoidal,2],AXIS["longitude",east,ORDER[1],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]],AXIS["latitude",north,ORDER[2],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]]' to `ENGCRS["Cartesian (Meter)",EDATUM["Unknown engineering datum"],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]' (GDAL error 6)
#> Warning: [extract] transforming vector data to the CRS of the raster
#> Warning: number of rows of result is not a multiple of vector length (arg 1)
#> Error in idx[, 2]: subscript out of bounds

extractRange(r, xy, first=1:2, last=3:2)
#> Warning: Cannot find coordinate operations from `GEOGCRS["unknown",DATUM["World Geodetic System 1984",ELLIPSOID["WGS 84",6378137,298.257223563,LENGTHUNIT["metre",1]],ID["EPSG",6326]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]],CS[ellipsoidal,2],AXIS["longitude",east,ORDER[1],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]],AXIS["latitude",north,ORDER[2],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]]' to `ENGCRS["Cartesian (Meter)",EDATUM["Unknown engineering datum"],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]' (GDAL error 6)
#> Warning: [extract] transforming vector data to the CRS of the raster
#> Warning: number of rows of result is not a multiple of vector length (arg 1)
#> Error in first[i]:last[i]: argument of length 0
extractRange(r, xy, first=1:2, last=3:2, lyr_fun=sum)
#> Warning: Cannot find coordinate operations from `GEOGCRS["unknown",DATUM["World Geodetic System 1984",ELLIPSOID["WGS 84",6378137,298.257223563,LENGTHUNIT["metre",1]],ID["EPSG",6326]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]],CS[ellipsoidal,2],AXIS["longitude",east,ORDER[1],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]],AXIS["latitude",north,ORDER[2],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]]' to `ENGCRS["Cartesian (Meter)",EDATUM["Unknown engineering datum"],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]' (GDAL error 6)
#> Warning: [extract] transforming vector data to the CRS of the raster
#> Warning: number of rows of result is not a multiple of vector length (arg 1)
#> Error in first[i]:last[i]: argument of length 0
```

# Get or set a coordinate reference system

Get or set the coordinate reference system (CRS), also referred to as a
"projection", of a SpatRaster or SpatVector.

Setting a new CRS does not change the data itself, it just changes the
label. So you should only set the CRS of a dataset (if it does not come
with one) to what it \*is\*, not to what you would \*like it to be\*.
See [`project`](https://rspatial.github.io/terra/reference/project.md)
to \*transform\* an object from one CRS to another.

## Usage

``` r
# S4 method for class 'SpatRaster'
crs(x, proj=FALSE, describe=FALSE, parse=FALSE)

# S4 method for class 'SpatVector'
crs(x, proj=FALSE, describe=FALSE, parse=FALSE)

# S4 method for class 'character'
crs(x, proj=FALSE, describe=FALSE, parse=FALSE)

# S4 method for class 'SpatRaster'
crs(x, warn = FALSE) <- value

# S4 method for class 'SpatVector'
crs(x, warn = FALSE) <- value
```

## Arguments

- x:

  SpatRaster or SpatVector

- proj:

  logical. If `TRUE` the crs is returned in PROJ-string notation

- describe:

  logical. If `TRUE` the name, EPSG code, and the name and extent of the
  area of use are returned if known

- warn:

  logical. If `TRUE`, a message is printed when the object already has a
  non-empty crs

- value:

  character string describing a coordinate reference system. This can be
  in a WKT format, as a \<authority:number\> code such as "EPSG:4326",
  or a PROJ-string format such as "+proj=utm +zone=12" (see Note)

- parse:

  logical. If `TRUE`, wkt parts are parsed into a vector (each line
  becomes an element)

## Note

Projections are handled by the PROJ/GDAL libraries. The PROJ developers
suggest to define a CRS with the WKT2 or \<authority\>:\<code\>
notation. It is not practical to define one's own custom CRS with WKT2,
and the the \<authority\>:\<code\> system only covers a handful of
(commonly used) CRSs. To work around this problem it is still possible
to use the deprecated PROJ-string notation (`+proj=...`) with one major
caveat: the datum should be WGS84 (or the equivalent NAD83) – if you
want to transform your data to a coordinate reference system with a
different datum. Thus as long as you use WGS84, or an ellipsoid instead
of a datum, you can safely use PROJ-strings to represent your CRS;
including to define your own custom CRS.

You can also set the crs to "local" to get an informal coordinate system
on an arbitrary Euclidean (Cartesian) plane with units in meter.

## Value

character or modified SpatRaster/Vector

## Examples

``` r
r <- rast()
crs(r)
#> [1] "GEOGCRS[\"WGS 84 (CRS84)\",\n    DATUM[\"World Geodetic System 1984\",\n        ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n            LENGTHUNIT[\"metre\",1]]],\n    PRIMEM[\"Greenwich\",0,\n        ANGLEUNIT[\"degree\",0.0174532925199433]],\n    CS[ellipsoidal,2],\n        AXIS[\"geodetic longitude (Lon)\",east,\n            ORDER[1],\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        AXIS[\"geodetic latitude (Lat)\",north,\n            ORDER[2],\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n    USAGE[\n        SCOPE[\"unknown\"],\n        AREA[\"World\"],\n        BBOX[-90,-180,90,180]],\n    ID[\"OGC\",\"CRS84\"]]"
crs(r, describe=TRUE, proj=TRUE)
#>             name authority  code  area             extent
#> 1 WGS 84 (CRS84)       OGC CRS84 World -180, 180, -90, 90
#>                                  proj
#> 1 +proj=longlat +datum=WGS84 +no_defs

crs(r) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84"
crs(r)
#> [1] "PROJCRS[\"unknown\",\n    BASEGEOGCRS[\"unknown\",\n        DATUM[\"Unknown based on WGS 84 ellipsoid\",\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1],\n                ID[\"EPSG\",7030]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8901]]],\n    CONVERSION[\"unknown\",\n        METHOD[\"Lambert Conic Conformal (2SP)\",\n            ID[\"EPSG\",9802]],\n        PARAMETER[\"Latitude of false origin\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8821]],\n        PARAMETER[\"Longitude of false origin\",-100,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8822]],\n        PARAMETER[\"Latitude of 1st standard parallel\",48,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8823]],\n        PARAMETER[\"Latitude of 2nd standard parallel\",33,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8824]],\n        PARAMETER[\"Easting at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8826]],\n        PARAMETER[\"Northing at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8827]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]]]"

# You can use epsg codes
crs(r)  <- "epsg:25831"
crs(r, describe=TRUE)$area
#> [1] "Europe between 0°E and 6°E: Andorra; Belgium - onshore and offshore; Denmark - offshore; Germany - offshore; Jan Mayen - offshore; Norway including Svalbard - onshore and offshore; Spain - onshore and offshore"

crs("epsg:25831", describe=TRUE)
#>                    name authority  code
#> 1 ETRS89 / UTM zone 31N      EPSG 25831
#>                                                                                                                                                                                                                area
#> 1 Europe between 0°E and 6°E: Andorra; Belgium - onshore and offshore; Denmark - offshore; Germany - offshore; Jan Mayen - offshore; Norway including Svalbard - onshore and offshore; Spain - onshore and offshore
#>                     extent
#> 1 0.00, 6.01, 37.00, 82.45
```

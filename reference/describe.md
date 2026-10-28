# describe

Describe the properties of spatial data in a file as generated with the
"GDALinfo" tool.

## Usage

``` r
# S4 method for class 'character'
describe(x, sds=FALSE, meta=FALSE, parse=FALSE, options="", print=FALSE, open_opt="")

# S4 method for class 'SpatRaster'
describe(x, source, ...)
```

## Arguments

- x:

  character. The name of a file with spatial data. Or a fully specified
  subdataset within a file such as `"NETCDF:\"AVHRR.nc\":NDVI"`

- sds:

  logical. If `TRUE` the description or metadata of the subdatasets is
  returned (if available)

- meta:

  logical. Get the file level metadata instead

- parse:

  logical. If `TRUE`, metadata for subdatasets is parsed into components
  (if `meta=TRUE`)

- options:

  character. A vector of valid options (if `meta=FALSE`) including
  "json", "mm", "stats", "hist", "nogcp", "nomd", "norat", "noct",
  "nofl", "checksum", "proj4", "listmdd", "mdd \<value\>" where
  \<value\> specifies a domain or 'all', "wkt_format \<value\>" where
  value is one of 'WKT1', 'WKT2', 'WKT2_2015', or 'WKT2_2018', "sd
  \<subdataset\>" where \<subdataset\> is the name or identifier of a
  sub-dataset. See <https://gdal.org/en/latest/programs/gdalinfo.html>.
  Ignored if `sds=TRUE`

- print:

  logical. If `TRUE`, print the results

- open_opt:

  character. Driver specific open options

- source:

  positive integer between 1 and `nsrc(x)`

- ...:

  additional arguments passed to the `describe<character>` method

## Value

character (invisibly, if `print=FALSE`)

## See also

[`ar_info`](https://rspatial.github.io/terra/reference/ar_info.md)

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra")
describe(f)
#>  [1] "Driver: GTiff/GeoTIFF"                                                    
#>  [2] "Files: /home/runner/work/_temp/Library/terra/ex/elev.tif"                 
#>  [3] "Size is 95, 90"                                                           
#>  [4] "Coordinate System is:"                                                    
#>  [5] "GEOGCRS[\"WGS 84\","                                                      
#>  [6] "    ENSEMBLE[\"World Geodetic System 1984 ensemble\","                    
#>  [7] "        MEMBER[\"World Geodetic System 1984 (Transit)\"],"                
#>  [8] "        MEMBER[\"World Geodetic System 1984 (G730)\"],"                   
#>  [9] "        MEMBER[\"World Geodetic System 1984 (G873)\"],"                   
#> [10] "        MEMBER[\"World Geodetic System 1984 (G1150)\"],"                  
#> [11] "        MEMBER[\"World Geodetic System 1984 (G1674)\"],"                  
#> [12] "        MEMBER[\"World Geodetic System 1984 (G1762)\"],"                  
#> [13] "        MEMBER[\"World Geodetic System 1984 (G2139)\"],"                  
#> [14] "        MEMBER[\"World Geodetic System 1984 (G2296)\"],"                  
#> [15] "        ELLIPSOID[\"WGS 84\",6378137,298.257223563,"                      
#> [16] "            LENGTHUNIT[\"metre\",1]],"                                    
#> [17] "        ENSEMBLEACCURACY[2.0]],"                                          
#> [18] "    PRIMEM[\"Greenwich\",0,"                                              
#> [19] "        ANGLEUNIT[\"degree\",0.0174532925199433]],"                       
#> [20] "    CS[ellipsoidal,2],"                                                   
#> [21] "        AXIS[\"geodetic latitude (Lat)\",north,"                          
#> [22] "            ORDER[1],"                                                    
#> [23] "            ANGLEUNIT[\"degree\",0.0174532925199433]],"                   
#> [24] "        AXIS[\"geodetic longitude (Lon)\",east,"                          
#> [25] "            ORDER[2],"                                                    
#> [26] "            ANGLEUNIT[\"degree\",0.0174532925199433]],"                   
#> [27] "    USAGE["                                                               
#> [28] "        SCOPE[\"Horizontal component of 3D system.\"],"                   
#> [29] "        AREA[\"World.\"],"                                                
#> [30] "        BBOX[-90,-180,90,180]],"                                          
#> [31] "    ID[\"EPSG\",4326]]"                                                   
#> [32] "Data axis to CRS axis mapping: 2,1"                                       
#> [33] "Origin = (5.741666666666666,50.191666666666663)"                          
#> [34] "Pixel Size = (0.008333333333333,-0.008333333333333)"                      
#> [35] "Metadata:"                                                                
#> [36] "  AREA_OR_POINT=Area"                                                     
#> [37] "Image Structure Metadata:"                                                
#> [38] "  COMPRESSION=LZW"                                                        
#> [39] "  INTERLEAVE=BAND"                                                        
#> [40] "Corner Coordinates:"                                                      
#> [41] "Upper Left  (   5.7416667,  50.1916667) (  5d44'30.00\"E, 50d11'30.00\"N)"
#> [42] "Lower Left  (   5.7416667,  49.4416667) (  5d44'30.00\"E, 49d26'30.00\"N)"
#> [43] "Upper Right (   6.5333333,  50.1916667) (  6d32' 0.00\"E, 50d11'30.00\"N)"
#> [44] "Lower Right (   6.5333333,  49.4416667) (  6d32' 0.00\"E, 49d26'30.00\"N)"
#> [45] "Center      (   6.1375000,  49.8166667) (  6d 8'15.00\"E, 49d49' 0.00\"N)"
#> [46] "Band 1 Block=95x43 Type=Int16, ColorInterp=Gray"                          
#> [47] "  Description = elevation"                                                
#> [48] "  Min=141.000 Max=547.000 "                                               
#> [49] "  Minimum=141.000, Maximum=547.000, Mean=-9999.000, StdDev=-9999.000"     
#> [50] "  NoData Value=-32768"                                                    
#> [51] "  Metadata:"                                                              
#> [52] "    STATISTICS_MAXIMUM=547"                                               
#> [53] "    STATISTICS_MEAN=-9999"                                                
#> [54] "    STATISTICS_MINIMUM=141"                                               
#> [55] "    STATISTICS_STDDEV=-9999"                                              
describe(f, meta=TRUE)
#> [1] "AREA_OR_POINT=Area"
#g <- describe(f, options=c("json", "nomd", "proj4"))
#head(g)
```

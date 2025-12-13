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
xy <- data.frame(x=c(50,80), y=c(30, 60))
extract(r, xy)
#>   ID red green blue
#> 1  1 149   158  215
#> 2  2  68    67   63
extract(r, xy, layer=c("red", "green"))
#>   ID value
#> 1  1   149
#> 2  2    67

extractRange(r, xy, first=1:2, last=3:2)
#> [[1]]
#>   red green blue
#> 1 149   158  215
#> 
#> [[2]]
#>   green
#> 2    67
#> 
extractRange(r, xy, first=1:2, last=3:2, lyr_fun=sum)
#> [1] 522  67
```

# Split a SpatRaster or SpatVector

Split a SpatRaster by layer, or a SpatVector by attributes. You can also
split the geometry of a SpatVector of polygon or lines with another
SpatVector of polygon or lines.

## Usage

``` r
# S4 method for class 'SpatRaster,ANY'
split(x, f)

# S4 method for class 'SpatVector,ANY'
split(x, f)

# S4 method for class 'SpatVector,SpatVector'
split(x, f, min_node_dist=10000)
```

## Arguments

- x:

  SpatRaster or SpatVector

- f:

  If `x` is a SpatRaster: a vector of the length `nlyr(x)`. If `x` is a
  SpatVector: one or more variable names, or a vector of the same length
  as `x`, or a list of such vectors. If `x` is a SpatVector of polygons,
  you can also use a SpatVector of lines or polygons to split the
  polygon geometries

- min_node_dist:

  postive number indicating the minimum node distance to use (in m) for
  longitude/latitude data. To ensure this minium distance between nodes,
  additional nodes are added as needed, to improve precision. See
  [`densify`](https://rspatial.github.io/terra/reference/densify.md)

## Value

list or SpatVector

## See also

[`segregate`](https://rspatial.github.io/terra/reference/segregate.md)

## Examples

``` r
## split layers
s <- rast(system.file("ex/logo.tif", package="terra"))   
y <- split(s, c(1,2,1))
sds(y)
#> class       : SpatRasterDataset 
#> subdatasets : 2 
#> dimensions  : 77, 101 (nrow, ncol)
#> nlyr        : 2, 1 
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source(s)   : logo.tif 

## split attributes
v <- vect(system.file("ex/lux.shp", package="terra"))
x <- split(v, "NAME_1")

## split geometries
v <- v[1:5,]
line <- vect(matrix(c(5.79, 6.22, 5.75, 6.1, 5.8, 
  50.14, 50.05, 49.88, 49.85, 49.71), ncol=2), "line")
s <- split(v, line)
```

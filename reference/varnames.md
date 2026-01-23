# variable and long variable names

Set or get names for each dataset (variable) in a SpatRasterDataset.

Each SpatRaster \_data source\_ can also have a variable name and a long
variable name. They are set when reading a file with possibly multiple
sub-datasets (e.g. netcdf or hdf5 format) into a single SpatRaster. Each
sub-dataset is a separate "data-source" in the SpatRaster. Note that
newly created or derived SpatRasters always have a single variable (data
source), and therefore the variable names are lost when processing a
multi-variable SpatRaster. Thus the variable names are mostly useful to
understand a SpatRaster created from some files and for managing
SpatRasterDatasets.

See `link{names}` for the more commonly used \_layer\_ names.

## Usage

``` r
# S4 method for class 'SpatRaster'
varnames(x)

# S4 method for class 'SpatRaster'
varnames(x) <- value

# S4 method for class 'SpatRaster'
longnames(x)

# S4 method for class 'SpatRaster'
longnames(x) <- value

# S4 method for class 'SpatRasterDataset'
varnames(x)

# S4 method for class 'SpatRasterDataset'
varnames(x) <- value

# S4 method for class 'SpatRasterDataset'
longnames(x)

# S4 method for class 'SpatRasterDataset'
longnames(x) <- value
```

## Arguments

- x:

  SpatRaster, SpatRasterDataset

- value:

  character (vector)

## Value

character

## Note

terra enforces neither unique nor valid names. See
[`make.unique`](https://rdrr.io/r/base/make.unique.html) to create
unique names and `{make.names}` to make syntactically valid names.

## Examples

``` r
s <- rast(ncols=5, nrows=5, nlyrs=3)
names(s) <- c("a", "b", "c")
x <- sds(s, s)
varnames(x) <- c("one", "two")
x
#> class       : SpatRasterDataset 
#> subdatasets : 2 
#> dimensions  : 5, 5 (nrow, ncol)
#> nlyr        : 3, 3 
#> resolution  : 72, 36  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory 
#> names       : one, two 
```

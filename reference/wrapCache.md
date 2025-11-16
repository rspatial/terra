# SpatRaster wrap with caching options

Use `wrap` to pack a SpatRaster with caching options. See
[`wrap`](https://rspatial.github.io/terra/reference/wrap.md) for the
general approach that is easier and better to use in most cases.

This method allows for specifying a folder, or filenames, to cache all
sources of a SpatRaster in a specific location (on disk).

## Usage

``` r
# S4 method for class 'SpatRaster'
wrapCache(x, filename=NULL, path=NULL, overwrite=FALSE, ...)
```

## Arguments

- x:

  SpatRaster

- filename:

  character. A single filename, or one filename per SpatRaster data
  source. If not `NULL`, the raster sources are saved in these files

- path:

  character. If not `NULL`, the path where raster sources will be saved.
  Ignored if filenames is not `NULL`

- overwrite:

  Should existing files be overwritten when `files` or `path` is not
  `NULL`? If this value is not `TRUE` or `FALSE`, only files that do not
  exist are created

- ...:

  Additional arguments for `writeRaster`. Only used for raster sources
  that are in memory, as other sources are cached by copying the files

## Value

PackedSpatRaster

## See also

[`wrap`](https://rspatial.github.io/terra/reference/wrap.md),
[`unwrap`](https://rspatial.github.io/terra/reference/wrap.md)

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)

x <- wrapCache(r, path=tempdir())
x
#> [1] "This is a PackedSpatRaster object. Use 'terra::unwrap()' to unpack it"
```

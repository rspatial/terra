# ar_info

Describe a multi-dimensional array (netcdf) file

## Usage

``` r
ar_info(x, what="describe", simplify=TRUE, filter=TRUE, array="")
```

## Arguments

- x:

  character. The name of a netcdf (or similar) raster file

- what:

  character that (partially) matches "describe", "arrays" or
  "dimensions"

- simplify:

  logical. If `TRUE` and `what="describe"`, simplify the output for
  readability

- filter:

  logical. If `TRUE` and `what="describe"` filter arrays that (probably)
  dimensions

- array:

  character. Required when `what="dimensions"`

## Value

character or data.frame

## See also

[`describe`](https://rspatial.github.io/terra/reference/describe.md)

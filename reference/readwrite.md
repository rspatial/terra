# Read from, or write to, file

Methods to read from or write chunks of values to or from a file. These
are low level methods for programmers. Use writeRaster if you want to
save an entire SpatRaster to file in one step. It is much easier to use.

To write chunks, begin by opening a file with `writeStart`, then write
values to it in chunks using the list that is returned by `writeStart`.
When writing is done, close the file with `writeStop`.

`blocks` only returns chunk size information. This can be useful when
reading, but not writing, raster data.

## Usage

``` r
# S4 method for class 'SpatRaster'
readStart(x)

# S4 method for class 'SpatRaster'
readStop(x)

# S4 method for class 'SpatRaster'
readValues(x, row=1, nrows=nrow(x), col=1, ncols=ncol(x), mat=FALSE, dataframe=FALSE, ...)

# S4 method for class 'SpatRaster,character'
writeStart(x, filename="", overwrite=FALSE, n=4, sources="", ...)

# S4 method for class 'SpatRaster'
writeStop(x)

# S4 method for class 'SpatRaster,vector'
writeValues(x, v, start, nrows)

# S4 method for class 'SpatRaster'
blocks(x, n=4)

fileBlocksize(x)
```

## Arguments

- x:

  SpatRaster

- filename:

  character. Output filename

- v:

  vector with cell values to be written

- start:

  integer. Row number (counting starts at 1) from where to start writing
  `v`

- row:

  positive integer. Row number to start from, should be between 1 and
  nrow(x)

- nrows:

  positive integer. How many rows?

- col:

  positive integer. Column number to start from, should be between 1 and
  ncol(x)

- ncols:

  positive integer. How many columns? Default is the number of columns
  left after the start column

- mat:

  logical. If `TRUE`, values are returned as a numeric matrix instead of
  as a vector, except when `dataframe=TRUE`. If any of the layers of `x`
  is a factor, the level index is returned, not the label. Use
  `dataframe=TRUE` to get the labels

- dataframe:

  logical. If `TRUE`, values are returned as a `data.frame` instead of
  as a vector (also if matrix is `TRUE`)

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- n:

  positive integer indicating how many copies the data may be in memory
  at any point in time. This is used to determine how many blocks
  (large) datasets need to be read

- sources:

  character. Filenames that may not be overwritten because they are used
  as input to the function. Can be obtained with `sources(x)`

- ...:

  For `writeStart`: additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

  For `readValues`: additional arguments for
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) (and thus only
  relevant when `dataframe=TRUE`)

## Value

`readValues` returns a vector, matrix, or data.frame

`writeStart` returns a list that can be used for processing the file in
chunks.

The other methods invisibly return a logical value indicating whether
they were successful or not. Their purpose is the side-effect of opening
or closing files.

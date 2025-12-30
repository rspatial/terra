# Temporary files

List and optionally remove temporary files created by the terra package.
These files are created when an output SpatRaster may be too large to
store in memory (RAM). This can happen when no filename is provided to a
function and when using functions where you cannot provide a filename.

Temporary files are automatically removed at the end of each R session
that ends normally. You can use `tmpFiles` to see the files in the
current sessions, including those that are orphaned (not connect to a
SpatRaster object any more) and from other (perhaps old) sessions, and
remove all the temporary files.

## Usage

``` r
tmpFiles(current=TRUE, orphan=FALSE, old=FALSE, remove=FALSE)
```

## Arguments

- current:

  logical. If `TRUE`, temporary files from the current R session are
  included

- orphan:

  logical. If `TRUE`, temporary files from the current R session that
  are no longer associated with a SpatRaster (if `current` is `TRUE`
  these are also included)

- old:

  logical. If `TRUE`, temporary files from other "R" sessions. Unless
  you are running multiple instances of R at the same time, these are
  from old (possibly crashed) R sessions and should be removed

- remove:

  logical. If `TRUE`, temporary files are removed

## Value

character

## See also

[`terraOptions`](https://rspatial.github.io/terra/reference/terraOptions.md)

## Examples

``` r
tmpFiles()
#> [1] "/tmp/RtmplSS70d/spat_22e576e9cdac_8933_DvggOkjkTOQeAAI.vrt"
```

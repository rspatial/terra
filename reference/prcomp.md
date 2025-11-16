# SpatRaster PCA with prcomp

Compute principal components for SpatRaster layers. This method may be
preferred to
[`princomp`](https://rspatial.github.io/terra/reference/princomp.md) for
its greater numerical accuracy. However, it is slower and for very large
rasters it can only be done with a sample. This may be good enough but
see [`princomp`](https://rspatial.github.io/terra/reference/princomp.md)
if you want to use all values. Unlike
[`princomp`](https://rspatial.github.io/terra/reference/princomp.md), in
this method the sample variances are used with `n-1` as the denominator.

## Usage

``` r
# S4 method for class 'SpatRaster'
prcomp(x, retx=TRUE, center=TRUE, scale.=FALSE, 
    tol=NULL, rank.=NULL, maxcell=Inf)
```

## Arguments

- x:

  SpatRaster

- retx:

  a logical value indicating whether the rotated variables should be
  returned

- center:

  a logical value indicating whether the variables should be shifted to
  be zero centered. Alternately, a vector of length equal the number of
  columns of x can be supplied. The value is passed to
  [`scale`](https://rspatial.github.io/terra/reference/scale.md)

- scale.:

  a logical value indicating whether the variables should be scaled to
  have unit variance before the analysis takes place. The default is
  FALSE for consistency with S, but in general scaling is advisable.
  Alternatively, a vector of length equal the number of columns of x can
  be supplied. The value is passed to
  [`scale`](https://rspatial.github.io/terra/reference/scale.md)

- tol:

  a value indicating the magnitude below which components should be
  omitted. (Components are omitted if their standard deviations are less
  than or equal to tol times the standard deviation of the first
  component.) With the default null setting, no components are omitted
  (unless `rank.` is specified less than `min(dim(x))`). Other settings
  for `tol` could be `tol = 0` or `tol = sqrt(.Machine$double.eps)`,
  which would omit essentially constant components

- rank.:

  optionally, a number specifying the maximal rank, i.e., maximal number
  of principal components to be used. Can be set as alternative or in
  addition to tol, useful notably when the desired rank is considerably
  smaller than the dimensions of the matrix

- maxcell:

  positive integer. The maximum number of cells to be used. If this is
  smaller than ncell(x), a regular sample of `x` is used

## Value

prcomp object

## Note

`prcomp` may change the layer names if they are not valid. See
[`make.names`](https://rdrr.io/r/base/make.names.html). In that case,
you will get a warning, and would need to also make the layer names of
`x` valid before using `predict`. Even better would be to change them
before calling `prcomp`.

## See also

[`princomp`](https://rspatial.github.io/terra/reference/princomp.md),
[`prcomp`](https://rdrr.io/r/stats/prcomp.html)

## Examples

``` r
f <- system.file("ex/logo.tif", package = "terra")
r <- rast(f)
pca <- prcomp(r)
x <- predict(r, pca)

# use "index" to get a subset of the components
p <- predict(r, pca, index=1:2)
```

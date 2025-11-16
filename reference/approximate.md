# Estimate values for cell values that are `NA` by interpolating between layers

approximate uses the `stats` function
[`approx`](https://rdrr.io/r/stats/approxfun.html) to estimate values
for cells that are `NA` by interpolation across layers. Layers are
considered equidistant, unless argument `z` is used, or `time(x)`
returns values that are not `NA`, in which case these values are used to
determine distance between layers.

For estimation based on neighboring cells see
[`focal`](https://rspatial.github.io/terra/reference/focal.md)

## Usage

``` r
# S4 method for class 'SpatRaster'
approximate(x, method="linear", yleft, yright,
            rule=1, f=0, ties=mean, z=NULL, NArule=1,filename="",  ...)
```

## Arguments

- x:

  SpatRaster

- method:

  specifies the interpolation method to be used. Choices are "linear" or
  "constant" (step function; see the example in
  [`approx`](https://rdrr.io/r/stats/approxfun.html)

- yleft:

  the value to be returned before a non-`NA` value is encountered. The
  default is defined by the value of rule given below

- yright:

  the value to be returned after the last non-`NA` value is encountered.
  The default is defined by the value of rule given below

- rule:

  an integer (of length 1 or 2) describing how interpolation is to take
  place at for the first and last cells (before or after any non-`NA`
  values are encountered). If rule is 1 then NAs are returned for such
  points and if it is 2, the value at the closest data extreme is used.
  Use, e.g., `rule = 2:1`, if the left and right side extrapolation
  should differ

- f:

  for method = "constant" a number between 0 and 1 inclusive, indicating
  a compromise between left- and right-continuous step functions. If y0
  and y1 are the values to the left and right of the point then the
  value is `y0*(1-f)+y1*f` so that `f = 0)` is right-continuous and
  `f = 1` is left-continuous

- ties:

  Handling of tied 'z' values. Either a function with a single vector
  argument returning a single number result or the string "ordered"

- z:

  numeric vector to indicate the distance between layers (e.g., depth).
  The default is `time(x)` if these are not `NA` or else `1:nlys(x)`

- NArule:

  single integer used to determine what to do when only a single layer
  with a non-`NA` value is encountered (and linear interpolation is not
  possible). The default value of 1 indicates that all layers will get
  this value for that cell; all other values do not change the cell
  values

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

` `[`focal`](https://rspatial.github.io/terra/reference/focal.md),
[`fillTime`](https://rspatial.github.io/terra/reference/fillTime.md)

## Examples

``` r
r <- rast(ncols=5, nrows=5)
r1 <- setValues(r, runif(ncell(r)))
r2 <- setValues(r, runif(ncell(r)))
r3 <- setValues(r, runif(ncell(r)))
r4 <- setValues(r, runif(ncell(r)))
r5 <- setValues(r, NA)
r6 <- setValues(r, runif(ncell(r)))
r1[6:10] <- NA
r2[5:15] <- NA
r3[8:25] <- NA
s <- c(r1,r2,r3,r4,r5,r6)
s[1:5] <- NA
x1 <- approximate(s)
x2 <- approximate(s, rule=2)
x3 <- approximate(s, rule=2, z=c(1,2,3,5,14,15))
```

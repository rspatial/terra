# Classify (or reclassify) cell values

Classify values of a SpatRaster. The function (re-)classifies groups of
values to other values.

The classification is done based on the argument `rcl`. You can classify
ranges by specifying a three-column matrix "from-to-becomes" or change
specific values by using a two-column matrix "is-becomes". You can also
supply a vector with "cuts" or the "number of cuts".

With "from-to-becomes" or "is-becomes" classification is done in the row
order of the matrix. Thus, if there are overlapping ranges or values,
the first time a number is within a range determines the
reclassification value.

With "cuts" the values are sorted, so that the order in which they are
provided does not matter.

## Usage

``` r
# S4 method for class 'SpatRaster'
classify(x, rcl, include.lowest=FALSE, right=TRUE, 
     others=NULL, brackets=TRUE, filename="", ...)
```

## Arguments

- x:

  SpatRaster

- rcl:

  matrix for classification. This matrix must have 1, 2 or 3 columns. If
  there are three columns, the first two columns are "from" "to" of the
  input values, and the third column "becomes" has the new value for
  that range.

  The two column matrix ("is", "becomes") can be useful for classifying
  integer values. In that case, the arguments `right` and
  `include.lowest` are ignored. For this case, an optimized hash-based
  lookup is used for improved performance. Floating point values are
  supported and are matched exactly.

  A single column matrix (or a vector) is interpreted as a set of cuts
  if there is more than one value. In that case the values are
  classified based on their location in-between the cut-values.

  If a single number is provided, that is used to make that number of
  cuts, at equal intervals between the lowest and highest values of the
  SpatRaster.

- include.lowest:

  logical, indicating if a value equal to the lowest value in `rcl` (or
  highest value in the second column, for `right=FALSE`) should be
  included.

- right:

  logical. If `TRUE`, the intervals are closed on the right (and open on
  the left). If `FALSE` they are open at the right and closed at the
  left. "open" means that the extreme value is \*not\* included in the
  interval. Thus, right-closed and left open is
  `(0,1] = {x | 0 < x <= 1}`. You can also close both sides with
  `right=NA`, that is only meaningful if you "from-to-becomes"
  classification with integers. For example to classify 1-5 -\> 1, 6-10
  -\> 2, 11-15 -\> 3. That may be easier to read and write than the
  equivalent 1-5 -\> 1, 5-10 -\> 2, 10-15 -\> 3 with `right=TRUE` and
  `include.lowest=TRUE`

- others:

  numeric. If not `NULL` all values that are not matched are set to this
  value. Otherwise they retain their original value.

- brackets:

  logical. If `TRUE`, intervals are have parenthesis or brackets around
  them to indicate whether they are open or closed. Only applies if
  `rcl` is a vector (or single column matrix)

- filename:

  character. Output filename

- ...:

  Additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## Author

Andrew Gene Brown, Robert J. Hijmans

## See also

[`subst`](https://rspatial.github.io/terra/reference/subst.md) for
simpler from-to replacement, and
[`clamp`](https://rspatial.github.io/terra/reference/clamp.md)

## Note

classify works with the "raw" values of categorical rasters, ignoring
the levels (labels, categories). To change the labels of categorical
rasters, use
[`subst`](https://rspatial.github.io/terra/reference/subst.md) instead.

For model-based classification see
[`predict`](https://rspatial.github.io/terra/reference/predict.md)

## Examples

``` r
r <- rast(ncols=10, nrows=10)
values(r) <- (0:99)/99

## from-to-becomes
# classify the values into three groups 
# all values >= 0 and <= 0.25 become 1, etc.
m <- c(0, 0.25, 1,
       0.25, 0.5, 2,
       0.5, 1, 3)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc1 <- classify(r, rclmat, include.lowest=TRUE)

## cuts
# equivalent to the above, but now a categorical SpatRaster is returned 
rc2 <- classify(r, c(0, 0.25, 0.5, 1), include.lowest=TRUE, brackets=TRUE)
freq(rc2)
#>   layer        value count
#> 1     1   [0 - 0.25]    25
#> 2     1 (0.25 - 0.5]    25
#> 3     1    (0.5 - 1]    50

## is-becomes 
x <- round(r*3)
unique(x)
#>   lyr.1
#> 1     0
#> 2     1
#> 3     2
#> 4     3
# replace 0 with NA
y <- classify(x, cbind(0, NA))
unique(y)
#>   lyr.1
#> 1     1
#> 2     2
#> 3     3

# multiple replacements
m <- rbind(c(2, 200), c(3, 300))
m
#>      [,1] [,2]
#> [1,]    2  200
#> [2,]    3  300

rcx1 <- classify(x, m)
unique(rcx1)
#>   lyr.1
#> 1     0
#> 2     1
#> 3   200
#> 4   300

rcx2 <- classify(x, m, others=NA)
unique(rcx2)
#>   lyr.1
#> 1   200
#> 2   300
```

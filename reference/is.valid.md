# Check or fix polygon or extent validity

Check the validity of polygons or attempt to fix it. Or check the
validity of a SpatExtent.

## Usage

``` r
# S4 method for class 'SpatVector'
is.valid(x, messages=FALSE, as.points=FALSE)

# S4 method for class 'SpatVector'
makeValid(x, buffer=FALSE)

# S4 method for class 'SpatExtent'
is.valid(x)
```

## Arguments

- x:

  SpatVector or SpatExtent

- messages:

  logical. If `TRUE` the error messages are returned

- as.points:

  logical. If `TRUE`, it is attempted to return locations where polygons
  are invalid as a SpatVector or points

- buffer:

  logical. If `TRUE` the zero-width buffer method is used to create
  valid polygons. Be careful when using this method because it may
  result in data loss. For example, only a single part of a
  self-intersecting may be preserved. See the example below

## Value

logical

## See also

[`topology`](https://rspatial.github.io/terra/reference/topology.md)

## Examples

``` r
w <- vect("POLYGON ((0 -5, 10 0, 10 -10, 0 -5))")
is.valid(w)
#> [1] TRUE

w <- vect("POLYGON ((0 -5, 10 0, 10 -10, 4 -2, 0 -5))")
is.valid(w)
#> [1] FALSE
is.valid(w, TRUE)
#>   valid                                                reason
#> 1 FALSE Self-intersection[4.54545454545454 -2.72727272727273]
v1 <- makeValid(w)
is.valid(v1)
#> [1] TRUE
v2 <- makeValid(w, buffer=TRUE)

plot(disagg(v1), col=c("light blue", "gray"))
# valid but incomplete
lines(v2, col="red", lwd=3)


plot(w)
points(cbind(4.54, -2.72), cex=2, col="red")


e <- ext(0, 1, 0, 1)
is.valid(e)
#> [1] TRUE

ee <- ext(0, 0, 0, 0)
is.valid(ee)
#> [1] TRUE
```

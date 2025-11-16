# Compute focal values with an iterating C++ function

Calculate focal values with a C++ function that iterates over cells to
speed up computations by avoiding an R loop (with apply).

See [`focal`](https://rspatial.github.io/terra/reference/focal.md) for
an easier to use method.

## Usage

``` r
# S4 method for class 'SpatRaster'
focalCpp(x, w=3, fun, ..., fillvalue=NA, 
    silent=TRUE, filename="", overwrite=FALSE, wopt=list())
```

## Arguments

- x:

  SpatRaster

- w:

  window. The window can be defined as one (for a square) or two numbers
  (row, col); or with an odd-sized weights matrix. See the Details
  section in
  [`focal`](https://rspatial.github.io/terra/reference/focal.md)

- fun:

  [`cppFunction`](https://rdrr.io/pkg/Rcpp/man/cppFunction.html) that
  iterates over cells. For C++ functions that operate on a single focal
  window, or for R functions use
  [`focal`](https://rspatial.github.io/terra/reference/focal.md)
  instead. The function must have at least three arguments. The first
  argument can have any name, but it must be a `Rcpp::NumericVector`,
  `Rcpp::IntegerVector` or a `std::vector<double>`. This is the
  container that receives the focal values. The other two arguments `ni`
  and `wi` must be of type `size_t`. `ni` represents the number of cells
  and `nw` represents the size of (number of elements in) the window

- ...:

  additional arguments to `fun`

- fillvalue:

  numeric. The value of the cells in the virtual rows and columns
  outside of the raster

- silent:

  logical. If `TRUE` error messages are printed that may occur when
  trying `fun` to determine the length of the returned value. This can
  be useful in debugging a `fun` that does not work

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- wopt:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`focal`](https://rspatial.github.io/terra/reference/focal.md),
[`focalValues`](https://rspatial.github.io/terra/reference/focalValues.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(Rcpp)
cppFunction( 
  'NumericVector sum_and_multiply(NumericVector x, double m, size_t ni, size_t nw) {
    NumericVector out(ni);
    // loop over cells
    size_t start = 0;
    for (size_t i=0; i<ni; i++) {
      size_t end = start + nw;
      // compute something for a window
      double v = 0;
      // loop over the values of a window
      for (size_t j=start; j<end; j++) {
        v += x[j];
      }
      out[i] = v * m;
      start = end;
    }
    return out;
  }'
)

nr <- nc <- 10
r <- rast(ncols=nc, nrows=nr, ext= c(0, nc, 0, nr))
values(r) <- 1:ncell(r)

raw <- focalCpp(r, w=3, fun=sum_and_multiply, fillvalue=0, m=10)

# same as
f1 <- focal(r, w=3, fun=sum, fillvalue=0) *10
all(values(f1) == values(raw))

# and as
ffun <- function(x, m) { sum(x) * m }
f2 <- focal(r, w=3, fun=ffun, fillvalue=0, m=10)


# You can also use an R function with focalCpp but this
# is not recommended 

R_sm_iter <- function(x, m, ni, nw) {
  out <- NULL
  for (i in 1:ni) {
    start <- (i-1) * nw + 1
    out[i] <- sum(x[start:(start+nw-1)]) * m
  }
  out
}

fr <- focalCpp(r, w=3, fun=R_sm_iter, fillvalue=0, m=10)

} # }
```

# Cell level regression

Run a regression model for each cell of a SpatRaster. The independent
variable can be defined either by:

- a `numeric` vector of length `nlyr(y)` – one covariate that varies
  across the layers but is constant across cells (e.g. time);

- a `SpatRaster` of the same number of layers as `y` – one covariate
  that varies both across layers and across cells;

- a `data.frame` with one row per layer of `y` and one column per
  predictor – multiple covariates that vary across layers but are
  constant across cells (e.g. an experimental design like NPK fertilizer
  levels).

## Usage

``` r
# S4 method for class 'SpatRaster,numeric'
regress(y, x, formula=y~x, na.rm=FALSE, cores=1, filename="", overwrite=FALSE, ...)

# S4 method for class 'SpatRaster,SpatRaster'
regress(y, x, formula=y~x, na.rm=FALSE, cores=1, filename="", overwrite=FALSE, ...)

# S4 method for class 'SpatRaster,data.frame'
regress(y, x, formula=NULL, na.rm=FALSE, cores=1, filename="", overwrite=FALSE, ...)
```

## Arguments

- y:

  SpatRaster. The dependent variable. Each layer is one observation per
  cell.

- x:

  The independent variable(s):

  - `SpatRaster` of the same number of layers as `y`;

  - `numeric` vector of length `nlyr(y)`;

  - `data.frame` with `nlyr(y)` rows and one column per predictor (no
    `NA`s allowed).

- formula:

  regression formula. For the `numeric` and `SpatRaster` methods the
  default is `y ~ x`; you can add additional terms such as `I(x^2)`. For
  the `data.frame` method the default (when `formula=NULL`) is `y ~ .`,
  expanding to all columns of `x`; you can also pass a custom formula
  whose right-hand side references the column names of `x` (the
  left-hand side is ignored).

- na.rm:

  logical. If `TRUE`, layers with `NA` response in a given cell are
  dropped before fitting. A cell whose remaining row count drops below
  the number of regression coefficients yields all `NA`s.

- cores:

  positive integer. If `cores > 1`, a 'parallel' package cluster with
  that many cores is created and used. You can also supply a cluster
  object.

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- ...:

  list with named options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster. One layer per regression coefficient, named after the
corresponding column of `model.matrix(formula, x)` (e.g. `(Intercept)`,
the names of the predictors, and any interaction or polynomial terms).

## Examples

``` r
s <- rast(system.file("ex/logo.tif", package="terra"))
x <- regress(s, 1:nlyr(s))

# data.frame method: a small NPK fertilizer experiment.
# Suppose `y` has one layer per (N, P, K) treatment.
if (FALSE) { # \dontrun{
NPK <- data.frame(
    N = c(0, 0, 0, 50, 50, 50, 100, 100, 100),
    P = c(0, 25, 50, 0, 25, 50, 0, 25, 50),
    K = c(0, 0, 0, 10, 10, 10, 20, 20, 20)
)
y <- rast(ncol=2, nrow=2, nlyr=nrow(NPK))
values(y) <- rep(1:nrow(NPK) * 200, each=4)

# Default formula y ~ N + P + K:
fit <- regress(y, NPK)

# Custom formula with interactions and a quadratic term:
fit2 <- regress(y, NPK, formula = yield ~ N + P + K + I(N^2) + N:P)
} # }
```

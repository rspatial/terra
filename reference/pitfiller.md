# Pit Filler

Fills pits (depressions with no outlet) in a elevation raster using the
PEM4PIT terrain erosion model.

## Usage

``` r
# S4 method for class 'SpatRaster'
pitfiller(x, pit, flowdir, niter = 10, lambda = 0,
                 deviation_type = "lad", max_iters=10^6, U = 1, D = 300, beta = 0.9,
                 theta_exp = 0.5, filename = "", ...)
```

## Arguments

- x:

  SpatRaster with elevation values

- pit:

  SpatRaster with pits identified with values \> 0. See
  [`pitfinder`](https://rspatial.github.io/terra/reference/pitfinder.md)

- flowdir:

  SpatRaster with flow direction. See
  [`flowDir`](https://rspatial.github.io/terra/reference/flowDir.md)

- niter:

  Number of iterations

- lambda:

  Deviation parameter

- deviation_type:

  Type of deviation. See
  [`flowDir`](https://rspatial.github.io/terra/reference/flowDir.md)

- max_iters:

  maximum iterations for drainage path starting points detection.See
  [`flowDir`](https://rspatial.github.io/terra/reference/flowDir.md)

- U:

  Uplift rate. See Details

- D:

  Diffusion coefficient. See Details

- beta:

  Erosion coefficient. See Details

- theta_exp:

  Exponent for contributing area. See Details

- filename:

  Character. Output filename

- ...:

  Additional arguments for writing files, as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

A
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.md)
with pits filled.

## Details

The PEM4PIT (Physical Erosion Model for PIT and flat areas correction)
equation \[Grimaldi at al, 2007\] corrects DEMs by simulating erosion
processes. It improves hydrologic modeling by enhancing flow direction
and contributing area representation. Function is still experimental.

The model is based on the steady-state sediment balance equation:

\$\$0 = U - \beta A^{\theta} \left( \frac{z - z_d}{\Delta l} \right) +
\frac{4D}{\Delta x^2} \left( \bar{z} - z \right)\$\$

Rearranging terms to isolate \\z\\:

\$\$\beta A^{\theta} \left( \frac{z - z_d}{\Delta l} \right) = U +
\frac{4D}{\Delta x^2} \left( \bar{z} - z \right)\$\$

Distributing terms:

\$\$\beta A^{\theta} \frac{z}{\Delta l} - \beta A^{\theta}
\frac{z_d}{\Delta l} = U + \frac{4D \bar{z}}{\Delta x^2} - \frac{4D
z}{\Delta x^2}\$\$

Combining like terms:

\$\$\beta A^{\theta} \frac{z}{\Delta l} + \frac{4D z}{\Delta x^2} = U +
\frac{4D \bar{z}}{\Delta x^2} + \beta A^{\theta} \frac{z_d}{\Delta
l}\$\$

Factoring out \\z\\:

\$\$z \left( \beta A^{\theta} \frac{1}{\Delta l} + \frac{4D}{\Delta x^2}
\right) = U + \frac{4D \bar{z}}{\Delta x^2} + \beta A^{\theta}
\frac{z_d}{\Delta l}\$\$

Solving for \\z\\:

\$\$z = \frac{U \Delta l \Delta x^2 + 4D \Delta l \bar{z} + \beta
A^{\theta} \Delta x^2 z_d}{\beta A^{\theta} \Delta x^2 + 4D \Delta
l}\$\$

## References

Grimaldi, S., Nardi, F., Di Benedetto, F., Istanbulluoglu, E., & Bras,
R. L. (2007). A physically-based method for removing pits in digital
elevation models. *Advances in Water Resources*, 30(10), 2151–2158.
[doi:10.1016/j.advwatres.2006.11.016](https://doi.org/10.1016/j.advwatres.2006.11.016)

Tucker, G. E., & Bras, R. L. (1998). Hillslope processes, drainage
density, and landscape morphology. *Water Resources Research*, 34(10),
2751–2764. [doi:10.1029/98WR01474](https://doi.org/10.1029/98WR01474)

Niemann, J. D., Bras, R. L., & Veneziano, D. (2003). A physically based
interpolation method for fluvially eroded topography. *Water Resources
Research*, 39, 1017.
[doi:10.1029/2001WR001050](https://doi.org/10.1029/2001WR001050)

Niemann, J. D., Bras, R. L., Veneziano, D., & Rinaldo, A. (2001).
Impacts of surface elevation on the growth and scaling properties of
simulated river networks. *Geomorphology*, 40(1–2), 37–55.
[doi:10.1016/S0169-555X(01)00036-8](https://doi.org/10.1016/S0169-555X%2801%2900036-8)

## Author

Emanuele Cordano

## See also

[`terrain`](https://rspatial.github.io/terra/reference/terrain.md),
[`watershed`](https://rspatial.github.io/terra/reference/watershed.md),
[`flowDir`](https://rspatial.github.io/terra/reference/flowDir.md),
[`pitfinder`](https://rspatial.github.io/terra/reference/pitfinder.md)

## Examples

``` r
elev <- rast(system.file("ex/elev.tif", package = "terra"))
elev <- project(elev, "epsg:32632")

lambda <- 0.5 ## try also 0 (default)
flowdir <- flowDir(elev, lambda = lambda)
pits <- pitfinder(flowdir, pits_on_boundary = FALSE)
elev2 <- pitfiller(x = elev, pit = pits,lambda=lambda)
#> Error in (function (cond) .Internal(C_tryCatchHelper(addr, 1L, cond)))(structure(list(message = "argument \"flowdir\" is missing, with no default",     call = .local(x, pit, ...)), class = c("evalError", "missingArgError", "error", "condition"))): error in evaluating the argument 'x' in selecting a method for function 'mask': argument "flowdir" is missing, with no default

flowdir2 <- terrain(elev2, "flowdir")
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'terrain': object 'elev2' not found
flowdir2 <- flowDir(elev2, lambda = lambda)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'flowDir': object 'elev2' not found
pits2 <- pitfinder(flowdir, pits_on_boundary = FALSE)
```

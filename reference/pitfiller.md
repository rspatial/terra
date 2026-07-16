# Pit Filler

Fills pits (depressions with no outlet) in a digital elevation model
(DEM) using the PEM4PIT terrain erosion model.

## Usage

``` r
# S4 method for class 'SpatRaster'
pitfiller(x, pit = NULL, flowdir = NULL, niter = 10, lambda = 0,
                 deviation_type = "lad", max_iters=10^6, U = 1, D = 300, beta = 0.9,
                 theta_exp = 0.5, filename = "", ...)
```

## Arguments

- x:

  SpatRaster with digital elevation model.

- pit:

  SpatRaster with pits. See
  [`pitfinder`](https://rspatial.github.io/terra/reference/pitfinder.md).

- flowdir:

  SpatRaster with flow direction or `NULL`. If `NULL`, it is calculated
  internally. See
  [`flowdirD8lad`](https://rspatial.github.io/terra/reference/flowdirD8ltd.md).

- niter:

  Number of iterations. Default is 10.

- lambda:

  Deviation parameter. Default is 0.

- deviation_type:

  Type of deviation. Default is `"lad"`. See
  [`flowdirD8lad`](https://rspatial.github.io/terra/reference/flowdirD8ltd.md).

- max_iters:

  maximum iterations for drainage path starting points detection.See
  [`flowdirD8lad`](https://rspatial.github.io/terra/reference/flowdirD8ltd.md).

- U:

  Uplift rate. See equation parameters \[Grimaldi at al, 2007\].

- D:

  Diffusion coefficient. See equation parameters \[Grimaldi at al,
  2007\].

- beta:

  Erosion coefficient. See equation parameters \[Grimaldi at al, 2007\].

- theta_exp:

  Exponent for contributing area. See equation parameters \[Grimaldi at
  al, 2007\].

- filename:

  Character. Output filename.

- ...:

  Additional arguments for writing files, as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md).

## Value

A
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.md)
object with pits filled.

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
[`flowdirD8lad`](https://rspatial.github.io/terra/reference/flowdirD8ltd.md),
[`pitfinder`](https://rspatial.github.io/terra/reference/pitfinder.md)

## Examples

``` r
#### PLANAR HILLSOPE WITH PIT EXAMPLE 

### PLANAR HILLSLOPE CREATION
# Parameters
res <- 100           # resolution in meters
size_km <- 1       # size of the area in kilometers
size_m <- size_km * 1000  # convert to meters

ncol <- size_m / res  # number of columns
nrow <- size_m / res  # number of rows

# Create an empty raster
r <- rast(nrows = nrow, ncols = ncol,
          xmin = 0, xmax = size_m, ymin = 0, ymax = size_m,
          resolution = res, crs = "")

# Slope angle in degrees
slope_deg <- 20
slope_rad <- slope_deg * pi / 180  # convert to radians

# Get Y coordinates of cell centers
y_coords <- yFromRow(r, 1:nrow)

# Compute elevation: elevationkg = distance * tan(slope)
elevation <- outer(y_coords, rep(1, ncol), function(y, x) y * tan(slope_rad))
# Assign elevation values to raster
values(r) <- as.vector(elevation)
elev0 <- r 

### PLANAR HILLSLOPE 
lambda=0##

flowdir0 <- flowdirD8lad(elev0, lambda = lambda)




plot(elev0,col=rev(terrain.colors(10)))
arrows_on_rast(flowdir0, unit="flowdir",col="black",code=2,length=0.1)

#> NULL

### PLANAR HILLSOPE WITH PIT 

elev1 <- elev0
elev1[47] <- elev1[47]-40 ##40
flowdir1 <- flowdirD8ltd(elev1, lambda = lambda)
flowdir1a <- terrain(elev1,"flowdir")
flowdir1a[flowdir1==0] <- 0
plot(elev1,col=rev(terrain.colors(10)))
##plot(elev1>as.numeric(elev1[47]))
##arrows_on_rast(flowdir1a, unit="flowdir",col="blue",code=2,length=0.1)
arrows_on_rast(flowdir1, unit="flowdir",col="black",code=2,length=0.1)

#> NULL

#### PIT DETECION AND FILLING 

elev <- elev1


flowdir <- flowdirD8lad(elev, lambda = lambda)
pits <- pitfinder(flowdir, pits_on_boundary = FALSE)
elev2 <- pitfiller(x = elev, pit = pits,lambda=lambda,niter=1000)
#> class       : SpatRaster
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 100, 100  (x, y)
#> extent      : 0, 1000, 0, 1000  (xmin, xmax, ymin, ymax)
#> coord. ref. : 
#> source(s)   : memory
#> name        : flowdir_lad_l=0
#> min value   :               0
#> max value   :               2


flowdir2 <- flowdirD8lad(elev2, lambda = lambda)
flowdir2a <- terrain(elev2, "flowdir")
flowdir2a[flowdir2==0] <- 0

pits2 <- pitfinder(flowdir2, pits_on_boundary = FALSE)

plot(elev2,col=rev(terrain.colors(10)))

arrows_on_rast(flowdir2a, unit="flowdir",col="blue",code=2,length=0.1)
#> NULL
arrows_on_rast(flowdir2, unit="flowdir",col="black",code=2,length=0.1)

#> NULL




#### LUX DIGITAL ELEVATION MODEL 
f <- system.file("ex/elev.tif", package = "terra")
elev <- rast(f) |> project(y = "epsg:32632")


lambda <- 0.5 ## try also 0 (default)
flowdir <- flowdirD8lad(elev, lambda = lambda)
pits <- pitfinder(flowdir, pits_on_boundary = FALSE)
elev2 <- pitfiller(x = elev, pit = pits,lambda=lambda)
#> class       : SpatRaster
#> size        : 111, 78, 1  (nrow, ncol, nlyr)
#> resolution  : 772.033, 772.033  (x, y)
#> extent      : 263811.2, 324029.8, 5479328, 5565024  (xmin, xmax, ymin, ymax)
#> coord. ref. : WGS 84 / UTM zone 32N (EPSG:32632)
#> source(s)   : memory
#> name        : flowdir_lad_l=0.5
#> min value   :                 0
#> max value   :               112
#> 
#> Exceeding number of iterations in d8ltd/d8lad flow directions computation
#> 
#> Exceeding number of iterations in d8ltd/d8lad flow directions computation
#> 
#> Exceeding number of iterations in d8ltd/d8lad flow directions computation
#> 
#> Exceeding number of iterations in d8ltd/d8lad flow directions computation
#> 
#> Exceeding number of iterations in d8ltd/d8lad flow directions computation
#> 
#> Exceeding number of iterations in d8ltd/d8lad flow directions computation
#> 
#> Exceeding number of iterations in d8ltd/d8lad flow directions computation
#> 
#> Exceeding number of iterations in d8ltd/d8lad flow directions computation
#> 
#> Exceeding number of iterations in d8ltd/d8lad flow directions computation
#> 
#> Exceeding number of iterations in d8ltd/d8lad flow directions computation
flowdir2 <- terrain(elev2, "flowdir")
flowdir2 <- flowdirD8lad(elev2, lambda = lambda)
pits2 <- pitfinder(flowdir, pits_on_boundary = FALSE)


```

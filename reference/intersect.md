# Intersection

You can intersect SpatVectors with each other or with a SpatExtent.
Intersecting points with points uses the extent of `y` to get the
intersection. Intersecting of points and lines is not supported because
of numerical inaccuracies with that. You can use
[`buffer`](https://rspatial.github.io/terra/reference/buffer.md), to
create polygons from lines and use these with intersect.

You can also intersect two SpatExtents.

When intersecting two SpatRasters these need to be aligned (have the
same origin and spatial resolution). The values of the returned
SpatRaster are `TRUE` where both input rasters have values, `FALSE`
where one has values, and `NA` in all other cells.

When intersecting a SpatExtent and a SpatRaster, the SpatExtent is first
aligned to the raster cell boundaries.

See [`crop`](https://rspatial.github.io/terra/reference/crop.md) for the
intersection of a SpatRaster with a SpatExtent (or the extent of a
SpatRaster or SpatVector) if you want a SpatRaster (not a SpatExtent) as
output.

See
[`is.related`](https://rspatial.github.io/terra/reference/relate.md)`(x, y, "intersects")`
to find out which geometries of a SpatVector intersect. You can
spatially subset a SpatVector with another one with
`x`[`[`](https://rspatial.github.io/terra/reference/subset_single.md)`y]`.

## Usage

``` r
# S4 method for class 'SpatVector,SpatVector'
intersect(x, y)

# S4 method for class 'SpatVector,SpatExtent'
intersect(x, y)

# S4 method for class 'SpatExtent,SpatVector'
intersect(x, y)

# S4 method for class 'SpatExtent,SpatExtent'
intersect(x, y)

# S4 method for class 'SpatRaster,SpatRaster'
intersect(x, y)

# S4 method for class 'SpatRaster,SpatExtent'
intersect(x, y)

# S4 method for class 'SpatExtent,SpatRaster'
intersect(x, y)
```

## Arguments

- x:

  SpatVector, SpatExtent, or SpatRaster

- y:

  SpatVector, SpatExtent, or SpatRaster

## Value

Same as `x`

## See also

[`union`](https://rspatial.github.io/terra/reference/union.md),
[`crop`](https://rspatial.github.io/terra/reference/crop.md),
[`relate`](https://rspatial.github.io/terra/reference/relate.md),
[`[`](https://rspatial.github.io/terra/reference/subset_single.md)

## Examples

``` r
e1 <- ext(-10, 10, -20, 20)
e2 <- ext(0, 20, -40, 5)
intersect(e1, e2)
#> SpatExtent : 0, 10, -20, 5 (xmin, xmax, ymin, ymax)

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
e <- ext(5.6, 6, 49.55, 49.7)
x <- intersect(v, e)

p <- vect(c("POLYGON ((5.8 49.8, 6 49.9, 6.15 49.8, 6 49.6, 5.8 49.8))", 
"POLYGON ((6.3 49.9, 6.2 49.7, 6.3 49.6, 6.5 49.8, 6.3 49.9))"), crs=crs(v))
values(p) <- data.frame(pid=1:2, area=expanse(p))

y <- intersect(v, p)

r <- s <- rast(ncol=5, nrow=5, xmin=1, xmax=5, ymin=1, ymax=5)
r[5:20] <- 5:20
s[11:20] <- 11:20
rs <- intersect(r, s)

u <- shift(r, .8)
us <- intersect(u, s)
```

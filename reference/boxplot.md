# Box plot of SpatRaster data

Box plot of layers in a SpatRaster

## Usage

``` r
# S4 method for class 'SpatRaster'
boxplot(x, y=NULL, maxcell=100000, ...)
```

## Arguments

- x:

  SpatRaster

- y:

  NULL or a SpatRaster. If `x` is a SpatRaster it used to group the
  values of `x` by "zone"

- maxcell:

  Integer. Number of cells to sample from datasets

- ...:

  additional arguments passed to
  `graphics::`[`boxplot`](https://rdrr.io/r/graphics/boxplot.html)

## Value

boxplot returns a list (invisibly) that can be used with
[`bxp`](https://rdrr.io/r/graphics/bxp.html)

## See also

[`pairs`](https://rspatial.github.io/terra/reference/pairs.md)`, `[`hist`](https://rspatial.github.io/terra/reference/hist.md)

## Examples

``` r
r1 <- r2 <- r3 <- rast(ncols=10, nrows=10)
set.seed(409)
values(r1) <- rnorm(ncell(r1), 100, 40)
values(r2) <- rnorm(ncell(r1), 80, 10)
values(r3) <- rnorm(ncell(r1), 120, 30)
s <- c(r1, r2, r3)
names(s) <- c("Apple", "Pear", "Cherry")

boxplot(s, notch=TRUE, col=c("red", "blue", "orange"), main="Box plot", ylab="random", las=1)


op <- par(no.readonly = TRUE)
par(mar=c(4,6,2,2))
boxplot(s, horizontal=TRUE, col="lightskyblue", axes=FALSE)
axis(1)
axis(2, at=0:3, labels=c("", names(s)), las=1, cex.axis=.9, lty=0)

par(op)

## boxplot with 2 layers
v <- vect(system.file("ex/lux.shp", package="terra"))
r <- rast(system.file("ex/elev.tif", package="terra"))
y <- rasterize(v, r, "NAME_2")
b <- boxplot(r, y)

bxp(b)
```

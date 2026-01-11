# Plot with leaflet

Plot a SpatRaster(Collection) or SpatVector(s) to make an interactive
leaflet map that is displayed in your browser.

The arguments of `plet` are similar to those of
[`plot`](https://rspatial.github.io/terra/reference/plot.md), making it
easier to use leaflet (if you also use plot).

## Usage

``` r
# S4 method for class 'SpatRaster'
plet(x, y=1, col, alpha=0.8, main=names(x), 
  tiles=c("Streets", "Esri.WorldImagery", "OpenTopoMap"), 
  wrap=TRUE, maxcell=500000, stretch=NULL, legend="bottomright", 
  shared=FALSE, panel=FALSE, collapse=TRUE, type=NULL, breaks=NULL,
  breakby="eqint", range=NULL, fill_range=FALSE, map=NULL, ...) 


# S4 method for class 'SpatRasterCollection'
plet(x, col, alpha=0.8, main=names(x), 
  tiles=c("Streets", "Esri.WorldImagery", "OpenTopoMap"), 
  wrap=TRUE, maxcell=500000, stretch=NULL, legend="bottomright", type=NULL, 
  breaks=NULL, breakby="eqint", range=NULL, fill_range=FALSE, map=NULL, ...)


# S4 method for class 'SpatVector'
plet(x, y="", col, main=y, cex=1, 
  lwd=2, lty=NULL, border="black", alpha=c(0.3, 1), popup=TRUE, label=FALSE,
  split=FALSE, tiles=c("Streets", "Esri.WorldImagery", "OpenTopoMap"), 
  wrap=TRUE, legend="bottomright", collapse=FALSE, type=NULL, breaks=NULL,
  breakby="eqint", sort=TRUE, reverse=FALSE, map=NULL, fill=NULL, ...)


# S4 method for class 'SpatVectorCollection'
plet(x, y="", col, main=y, cex=1, 
  lwd=2, lty=NULL, border="black", alpha=c(0.3, 1), popup=TRUE, label=FALSE, 
  tiles=c("Streets", "Esri.WorldImagery", "OpenTopoMap"), 
  wrap=TRUE, legend="bottomright", collapse=FALSE, type=NULL, breaks=NULL, 
  breakby="eqint", sort=TRUE, reverse=FALSE, map=NULL, fill=NULL, ...)


# S4 method for class 'leaflet'
lines(x, y, col, lwd=2, lty=NULL, alpha=1, ...)

# S4 method for class 'leaflet'
points(x, y, col, border=col, cex=1, lwd=2, lty=NULL, 
  alpha=c(.3, 1), label=1:nrow(y), popup=FALSE, ...)

# S4 method for class 'leaflet'
polys(x, y, col, lwd=2, lty=NULL, 
  border="black", alpha=c(0.3, 1), popup=TRUE, label=FALSE, fill=NULL, ...)
```

## Arguments

- x:

  SpatRaster, SpatVector, or leaflet object

- y:

  missing, or positive integer, or character (variable or layer name)
  indicating the layer(s) to be plotted. If `x` is a SpatRaster, you can
  select multiple layers

- col:

  character. Vector of colors or a color generating function. If `x` is
  a SpatVectorCollection, you can provide a list with colors and/or
  functions, with one list element for each SpatVector

- alpha:

  one or two numbers between 0 and 1 to set the transparency for lines
  (0 is transparent, 1 is opaque). The first number is to fill the
  points/lines/polygons, the second for the outline

- tiles:

  character or NULL. Names of background tile providers

- wrap:

  logical. if `TRUE`, tiles wrap around

- maxcell:

  positive integer. Maximum number of cells to use for the plot

- stretch:

  NULL or character ("lin" or "hist") to stretch RGB rasters. See
  [`plotRGB`](https://rspatial.github.io/terra/reference/plotRGB.md)

- legend:

  character to indicate the legend position ("bottomleft",
  "bottomright", "topleft" or "topright") or NULL to suppress the legend

- main:

  character. Title for the legend. The length should be 1 if `x` is a
  SpatVector and length nlyr(x) if `x` is a SpatVector

- shared:

  logical. Should the legend be the same for all rasters (if multiple
  layers of SpatRaster `x` are mapped)

- map:

  leaflet object

- ...:

  additional arguments for drawing points, lines, or polygons passed on
  the the relevant leaflet function

- border:

  character. Color for the polygon borders

- collapse:

  logical. Should the layers "control" panel be collapsed?

- split:

  logical. If `TRUE` a check-box is created to toggle each value in `y`
  (If `x` is a SpatVector)

- cex:

  numeric. point size magnifier. See
  [`par`](https://rdrr.io/r/graphics/par.html)

- lwd:

  numeric. line-width. See [`par`](https://rdrr.io/r/graphics/par.html)

- lty:

  character to specify a "dash-array". For example "3 5" indicates three
  pixels lines with five pixel gaps

- popup:

  logical. Should pop-ups be created?

- label:

  logical. Should mouse-over labels be added?

- panel:

  logical. Should SpatRaster layers be shown as a panel"

- type:

  character. Type of map/legend. One of "classes", or "interval". If not
  specified, the type is chosen based on the data. Use "" to suppress
  the legend

- breaks:

  numeric. Either a single number to indicate the number of breaks
  desired, or the actual breaks. When providing this argument, the
  default legend becomes "interval"

- breakby:

  character or function. Either "eqint" for equal interval breaks,
  "cases" for equal quantile breaks. If a function is supplied it should
  take a single argument (a vector of values) and create groups

- sort:

  logical. If `TRUE` legends with character values are sorted. You can
  also supply a vector of the unique values, in the order in which you
  want them to appear in the legend

- range:

  numeric. minimum and maximum values to be used for the continuous
  legend. You can use `NA` for one of these to only set the minimum or
  maximum value

- fill_range:

  logical. If `TRUE`, values outside of `range` get the colors of the
  extreme values; otherwise they get colored as `NA`

- reverse:

  logical. If `TRUE`, the legends order is reversed

- fill:

  do not use. Will be removed

## See also

[`plot`](https://rspatial.github.io/terra/reference/plot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (require(leaflet) && (packageVersion("leaflet") > "2.1.1")) {

v <- vect(system.file("ex/lux.shp", package="terra"))
p <- spatSample(as.polygons(v, ext=T), 30)
values(p) = data.frame(id=11:40, name=letters[1:30])

m <- plet(v, "NAME_1", tiles="", border="blue")
m <- points(m, p, col="red", cex=2, popup=T)
lines(m, v, lwd=1, col="white")

plet(v, "NAME_1", split=TRUE, alpha=.2) |> 
  points(p, col="white", border="red", cex=12, popup=TRUE, lwd=3, lty="1 4",
     clusterOptions = leaflet::markerClusterOptions())

s <- svc(v, p)
names(s) <- c("the polys", "set of points")
plet(s, col=c("red", "blue"), lwd=1)


r <- rast(system.file("ex/elev.tif", package="terra"))
plet(r, main="Hi\nthere", tiles=NULL) |> lines(v, lwd=1)

plet(r, tiles="OpenTopoMap") |> lines(v, lwd=2, col="blue")

x <- c(r, 50*classify(r, 5))
names(x) <- c("first", "second")

# each their own legend
plet(x, 1:2, collapse=FALSE) |> lines(v, lwd=2, col="blue", lty="5,5")

# shared legend
plet(x, 1:2, shared=TRUE, collapse=FALSE) |> lines(v, lwd=2, col="blue")

}} # }
```

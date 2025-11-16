# Change values in-place

These "in-place" replacement methods assign new value to an object
without making a copy. That is efficient, but if there is a copy of the
object that you made by standard assignment (e.g. with `y <- x`), that
copy is also changed.

`set.names` is the in-place replacement version of `names<-`.

`set.ext` is the in-place replacement version of `ext<-`

`set.values` is the in-place replacement version of `[<-`.

`set.cats` is the in-place replacement version of
[`categories`](https://rspatial.github.io/terra/reference/factors.md)

`set.crs` is the in-place replacement version of `crs<-`

`set.window` is the in-place replacement version of `window<-`

## Usage

``` r
# S4 method for class 'SpatRaster'
set.names(x, value, index=1:nlyr(x), validate=FALSE)
# S4 method for class 'SpatRasterDataset'
set.names(x, value, index=1:length(x), validate=FALSE)
# S4 method for class 'SpatVector'
set.names(x, value, index=1:ncol(x), validate=FALSE)

# S4 method for class 'SpatRaster'
set.ext(x, value)
# S4 method for class 'SpatVector'
set.ext(x, value)

# S4 method for class 'SpatRaster'
set.crs(x, value)
# S4 method for class 'SpatVector'
set.crs(x, value)

# S4 method for class 'SpatRaster'
set.values(x, cells, values, layer=0)
# S4 method for class 'SpatRasterDataset'
set.values(x)

# S4 method for class 'SpatRaster'
set.cats(x, layer=1, value, active=1)

# S4 method for class 'SpatRaster'
set.RGB(x, value, type="rgb")
```

## Arguments

- x:

  SpatRaster

- value:

  character for `set.names`. For `set.cats`: a data.frame with columns
  (value, category) or vector with category names. For `set.RGB` 3 or 4
  numbers indicating the RGB(A) layers

- index:

  positive integer indicating layer(s) to assign a name to

- validate:

  logical. Make names valid and/or unique?

- cells:

  cell numbers or missing

- values:

  replacement values or missing to load all values into memory

- layer:

  positive integer(s) indicating to which layer(s) to you want to assign
  these categories or to which you want to set these values. A number \<
  1 indicates "all layers"

- active:

  positive integer indicating the active category (column number in
  `value`, but not counting the first column

- type:

  character. The color space. One of "rgb" "hsv", "hsi" and "hsl"

## Value

logical (invisibly)

## Examples

``` r
s <- rast(ncols=5, nrows=5, nlyrs=3)
x <- s
names(s)
#> [1] "lyr.1" "lyr.2" "lyr.3"
names(s) <- c("a", "b", "c")
names(s)
#> [1] "a" "b" "c"
names(x)
#> [1] "lyr.1" "lyr.2" "lyr.3"

x <- s
set.names(s, c("e", "f", "g"))
names(s)
#> [1] "e" "f" "g"
names(x)
#> [1] "e" "f" "g"

set.ext(x, c(0,180,0,90))

f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)

#values from file to memory
set.values(r)
#> class       : SpatRaster 
#> size        : 90, 95, 1  (nrow, ncol, nlyr)
#> resolution  : 0.008333333, 0.008333333  (x, y)
#> extent      : 5.741667, 6.533333, 49.44167, 50.19167  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#> source(s)   : memory
#> varname     : elev 
#> name        : elevation 
#> min value   :       141 
#> max value   :       547 

# change values
set.values(r, 1:1000, 900)
```

# Plot a SpatExtent

Plot a SpatExtent. Use
[`lines`](https://rspatial.github.io/terra/reference/lines.md) to add a
SpatExtent to an existing map.

See [`plot`](https://rspatial.github.io/terra/reference/plot.md) for
plotting other object types.

## Usage

``` r
# S4 method for class 'SpatExtent,missing'
plot(x, y, ...)
```

## Arguments

- x:

  SpatExtent

- y:

  missing

- ...:

  additional graphical arguments for lines

## See also

[`plot`](https://rspatial.github.io/terra/reference/plot.md)

## Examples

``` r
r <- rast()
plot(ext(r))
```

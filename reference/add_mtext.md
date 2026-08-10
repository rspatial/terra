# add a margin text

Similar to [`mtext`](https://rdrr.io/r/graphics/mtext.html) allowing
adding a text to the margins of a map. This function uses the margins
around the mapped area; not the margins that R would use.

## Usage

``` r
add_mtext(text, side=3, line=0, adj=0.5, ...)
```

## Arguments

- text:

  character or expression vector specifying the text to be written

- side:

  integer indicating the margin to use (1=bottom, 2=left, 3=top,
  4=right)

- line:

  numeric to move the text in or outwards

- adj:

  numeric between 0 and 1 to align the text along the margin. Use 0 for
  left (sides 1 and 3) or bottom (sides 2 and 4), 1 for right or top,
  and 0.5, the default, to center it

- ...:

  arguments passed to
  [`text`](https://rspatial.github.io/terra/reference/text.md)

## See also

[`add_legend`](https://rspatial.github.io/terra/reference/legend.md),
[`add_grid`](https://rspatial.github.io/terra/reference/grid.md),
[`add_box`](https://rspatial.github.io/terra/reference/box.md)

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)

plot(r, axes=FALSE, legend=FALSE)
add_box()
for (i in 1:4) add_mtext("margin text", i, cex=i, col=i, line=2-i)


# align the text along the margin
plot(r, axes=FALSE, legend=FALSE)
add_box()
for (i in 1:4) {
  add_mtext("start", i, adj=0, col="blue")
  add_mtext("middle", i, adj=0.5, cex=1.5)
  add_mtext("end", i, adj=1, col="red")
}

```

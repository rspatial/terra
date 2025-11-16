# bestMatch

Determine for each grid cell which reference it is most similar to. A
reference consists of a SpatVector with reference locations, or a
data.frame or matrix in which each column matches a layer name in the
SpatRaster.

Similarity is computed with the mean absolute or the mean squared
differences between the cell and the reference, or with an alternative
function you provide. It may be important to first scale the input.

## Usage

``` r
# S4 method for class 'SpatRaster,SpatVector'
bestMatch(x, y, labels=NULL, fun="squared", ..., 
    filename="", overwrite=FALSE, wopt=list())

# S4 method for class 'SpatRaster,data.frame'
bestMatch(x, y, labels=NULL, fun="squared", ..., 
    filename="", overwrite=FALSE, wopt=list())

# S4 method for class 'SpatRaster,matrix'
bestMatch(x, y, labels=NULL, fun="squared", ..., 
    filename="", overwrite=FALSE, wopt=list())
```

## Arguments

- x:

  SpatRaster

- y:

  SpatVector, data.frame or matrix

- labels:

  character. labels that correspond to each class (row in `y`

- fun:

  character. One of "abs" for the mean absolute difference, or "squared"
  for the mean squared difference. Or a true function like
  terra:::match_sqr

- ...:

  additional arguments passed to `fun`. For the built-in functions this
  can be `na.rm=TRUE`

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- wopt:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## Examples

``` r
f <- system.file("ex/logo.tif", package = "terra")
r <- rast(f)

# locations of interest 
pts <- vect(cbind(c(25.25, 34.324, 43.003), c(54.577, 46.489, 30.905)))
pts$code <- LETTERS[1:3]

plot(r)
points(pts, pch=20, cex=2, col="red")
text(pts, "code", pos=4, halo=TRUE)


x <- scale(r)

s1 <- bestMatch(x, pts, labels=pts$code)
plot(s1)


# same result
e <- extract(x, pts, ID=FALSE)
s2 <- bestMatch(x, e, labels=c("Ap", "Nt", "Ms"))
```

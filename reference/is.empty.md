# Check if a SpatExtent or SpatVector is empty

An empty SpatExtent has no area

An empty SpatVector has no geometries.

## Usage

``` r
# S4 method for class 'SpatExtent'
is.empty(x)

# S4 method for class 'SpatVector'
is.empty(x)
```

## Arguments

- x:

  SpatVector or SpatExtent

## Value

logical

## Examples

``` r
e <- ext(0,0,0,0)
is.valid(e)
#> [1] TRUE
is.empty(e)
#> [1] TRUE

v <- vect()
is.valid(v)
#> logical(0)
is.empty(v)
#> [1] TRUE
```

# Memory available and needed

`mem_info` prints the amount of RAM that is required and available to
process a SpatRaster.

`free_RAM` returns the amount of RAM that is available

## Usage

``` r
mem_info(x, n=1, print=TRUE)

free_RAM()
```

## Arguments

- x:

  SpatRaster

- n:

  positive integer. The number of copies of `x` that are needed

- print:

  logical. print memory info?

## Value

free_RAM returns the amount of available RAM in kilobytes

## Examples

``` r
mem_info(rast())
#> 
#> ------------------------
#> Memory (GB) 
#> ------------------------
#> check threshold : 1 (memmin)
#> available       : 14.35
#> allowed (50%)   : 7.17
#> needed (n=1)    : 0
#> ------------------------
#> proc in memory  : TRUE
#> nr chunks       : 1
#> ------------------------

free_RAM()
#> [1] 15043032
```

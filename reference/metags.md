# Set or get metadata

You can set arbitrary metadata to (layers of) a SpatRaster using
"name=value", or "domain:name=value" tags or a two (name, value) or
three column (name, value, domain) matrix or data.frame.

## Usage

``` r
# S4 method for class 'SpatRaster'
metags(x, layer = NULL, domain = "") <- value

# S4 method for class 'SpatRaster'
metags(x, layer=NULL, name=NULL)

# S4 method for class 'SpatRasterDataset'
metags(x, dataset = NULL) <- value

# S4 method for class 'SpatRasterDataset'
metags(x, dataset=NULL, name=NULL)
```

## Arguments

- x:

  SpatRaster

- layer:

  NULL, positive integer or character. If the value is NULL, the tags
  assigned or returned are for the SpatRaster. Otherwise for the layer
  number(s) or name(s)

- domain:

  character. Only used if not specified by `value`. Use "" for the
  default domain. Depending on the file format used this may the only
  domain supported when writing files

- name:

  character

- value:

  character of "name=value" or two-column (name, value) or three-column
  (name, value, domain) matrix or data.frame

- dataset:

  NULL, positive integer or character. If the value is NULL, the tags
  assigned or returned are for the
  SpatRasterDataset/SpatRasterCollection. Otherwise for the datset
  number(s) or name(s)

## Value

SpatRaster (`metags<-`), or data.frame

## Examples

``` r
r <- rast(ncol=5, nrow=5)
m <- cbind(c("one", "two", "three"), c("ABC", "123", "hello"))
metags(r) <- m
metags(r)
#>    name value domain
#> 1   one   ABC       
#> 2   two   123       
#> 3 three hello       

metags(r) <- c("another_tag=another_value", "one more=this value")
metags(r)
#>          name         value domain
#> 1         one           ABC       
#> 2         two           123       
#> 3       three         hello       
#> 4 another_tag another_value       
#> 5    one more    this value       

metags(r) <- cbind("test", "this", "mydomain")
metags(r)
#>          name         value   domain
#> 1         one           ABC         
#> 2         two           123         
#> 3       three         hello         
#> 4 another_tag another_value         
#> 5    one more    this value         
#> 6        test          this mydomain

metags(r, name="two")
#>   name value domain
#> 2  two   123       

# remove a tag
metags(r) <- cbind("one", "")
metags(r) <- "two="
metags(r)
#>          name         value   domain
#> 1       three         hello         
#> 2 another_tag another_value         
#> 3    one more    this value         
#> 4        test          this mydomain

# remove all tags
metags(r) <- NULL
metags(r)
#> [1] name   value  domain
#> <0 rows> (or 0-length row.names)
```

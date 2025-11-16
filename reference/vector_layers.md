# List or remove layers from a vector file

List or remove layers from a vector file that supports layers such as
GPGK

## Usage

``` r
vector_layers(filename, delete="", return_error=FALSE)
```

## Arguments

- filename:

  character. filename

- delete:

  character. layers to be deleted (ignored if the value is `""`

- return_error:

  logical. If `TRUE`, an error occurs if some layers cannot be deleted.
  Otherwise a warning is given


# version 1.2-x

Under development

## new

- `trim` has a new argument `value` that allows trimming rows and columns with other values than the default `NA`.
- `rapp` has a new argument `clamp` that allows clamping start and end values to 1-nlyr(x), avoiding that all values are considered `NA`.

## bug fixes

- `extract` did not work for horizontal or vertical lines as their extent was considered invalid (bug reported by Monika Tomaszewska)
- `autocor` did not handle NA values [#192](https://github.com/rspatial/terra/issues/192)
- `nearest` now works for angular coordinates


# version 1.1-17

## major changes 

* `c-SpatVector-method` now returns a list. `rbind` is used to append SpatVector objects.
* overhaul of handling of factors. `rats` has been removed, and `levels` and `cats` have changed.


# version 1.1-4

- No news for this or earlier versions

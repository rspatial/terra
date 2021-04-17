
# version 1.2-x


## new

- `trim` has a new argument `value` that allows trimming rows and columns with other values than the default `NA`.
- `rapp` has a new argument `clamp` that allows clamping start and end values to `1:nlyr(x)`, avoiding that all values are considered `NA`.
- `spatSample,SpatRaster-method` has new arguments `as.points` and `values`. Getting values, cells and coordinates is no longer mutually exclusive. In response to [#191](https://github.com/rspatial/terra/issues/191). Requested by Agustin Lobo


## bug fixes

- `extract` did not work for horizontal or vertical lines as their extent was considered invalid. Reported by Monika Tomaszewska.
- `autocor` did not handle NA values [#192](https://github.com/rspatial/terra/issues/192). Reported by Laurence Hawker.
- `nearest` now works for angular coordinates
- The unit of `slope` in `terrain` was not correct (the tangent was returned instead of the slope) [#196](https://github.com/rspatial/terra/issues/196). Reported by Sven Alder.


# version 1.1-17

## major changes 

* `c-SpatVector-method` now returns a list. `rbind` is used to append SpatVector objects.
* overhaul of handling of factors. `rats` has been removed, and `levels` and `cats` have changed.


# version 1.1-4

- No news for this or earlier versions

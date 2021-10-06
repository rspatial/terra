# version 1.4-9

## name changes

To avoid name conflicts with `sp` (via `raster`) `disaggregate` is now called `disagg` and `bbox,SpatRaster` and `bbox,SpatVector` have been removed (but could be resurrected in `raster` or under another name).

## enhancements

- `project` and `resample` now choose the resampling method based on the first layer, using "near" for categorical data. Thanks to Matthew Lewis [#355](https://github.com/rspatial/terra/pull/355)

## bug fixes

- `hist` failed with small samples. Issue [#356](https://github.com/rspatial/terra/issues/356) by Martin Queinnec


# version 1.4-7

## note

`terra` no longer depends on `raster`. To avoid name clashes between these two packages, and to allow replacing methods from `rgeos` and `rgdal` in `raster`, a future version of `raster` will depend on `terra` instead. 


## enhancements

- `freq` has a new argument `usenames`. See issue [#309](https://github.com/rspatial/terra/issues/309) by Bappa Das
- `rast,character-method` has a new argument `opts` that can be used to pass GDAL open options. See issue [#314](https://github.com/rspatial/terra/issues/314)
- `rast,SpatRaster-method` now takes arguments `names` and `vals`. See issue [#323](https://github.com/rspatial/terra/issues/323) by Dongdong Kong
- `crs<-` now warns if an unsupported datum is used. See issue [#317](https://github.com/rspatial/terra/issues/317)
- `spatSample` now returns factor values if a SpatRaster layer is.factor except when using `as.df=FALSE`
- new method `origin<-` to set the origin of a SpatRaster. See issue [#326](https://github.com/rspatial/terra/issues/326) by Jakub Nowosad
- `crs` has a new argument `parse`. See [#344](https://github.com/rspatial/terra/issues/344) 
- `plot-SpatRaster,missing` has a new argument `reset=FALSE` that allows resetting the par()$mar parameters after plotting. See issue [#340](https://github.com/rspatial/terra/issues/340) by Derek Friend
- `crds` has a new argument `na.rm`. See [#338](https://github.com/rspatial/terra/issues/338) by Kodi Arfer 
- `show(Spat*)` now prints the name and EPSG code of a crs if available. See [#317](https://github.com/rspatial/terra/issues/317) by Jakub Nowosad


## bug fixes 

- `plotRGB` failed if there were `NA`s. Issue [#308](https://github.com/rspatial/terra/issues/308) by Jakub Nowosad
- `writeVector` crashed R when used with a SpatVector with no geometries. Reported by Timothy White in issue [#319](https://github.com/rspatial/terra/issues/319)
- `summary,SpatRaster-method` now returns counts for the classes (instead of a numerical summary of the indices) [#324](https://github.com/rspatial/terra/issues/324) by Jakub Nowosad
- `tapp` with a character index now returns a SpatRaster with the correct names [#345](https://github.com/rspatial/terra/issues/345) by Stuart Brown 
- `rasterize` with a character variable now adds the ID column to the categories [#337](https://github.com/rspatial/terra/issues/337) by Tate Brasel
- `cellSize` now masks values in all cases (when requested with `mask=TRUE`). Issue [#339](https://github.com/rspatial/terra/issues/339) by Jean-Luc Dupouey
- `buffer,SpatVector` no longer treats lines like polygons [#332](https://github.com/rspatial/terra/issues/332) by Márcia Barbosa
- `plot` now passes the layer index to `fun` [#310](https://github.com/rspatial/terra/issues/310) by Ben Tupper
- the `to_id` in `nearest` was sometimes wrong. See [#328](https://github.com/rspatial/terra/issues/328) by Shawn Ligocki
- better support for ESRI value attribute tables (VAT). See this [SO question]( https://stackoverflow.com/q/69385928/635245)
- `focal` did not reset initial values for NA cells when processing chunks. [#312](https://github.com/rspatial/terra/issues/312) by Jeffrey Evans
- `focal` could run out of memory when using a large window and user-defined function, and was inexact at the chunk boundary [#347](https://github.com/rspatial/terra/issues/347)
- `zonal` with `as.raster=TRUE` failed for categorical SpatRasters [#348](https://github.com/rspatial/terra/issues/348) by Jakub Nowosad



# version 1.3-22

## enhancements

- if `time(x) <- d` is set with a `Date` class object, `time(x)` now returns a `Date` object instead of a `POSIXct` object. Issue [#256](https://github.com/rspatial/terra/issues/256) raised by Mauricio Zambrano-Bigiarini
- The UTF-8 encoding of character attributes of a SpatVector is now declared such that they display correctly in R. See issue [#258](https://github.com/rspatial/terra/issues/258) by AGeographer. Also implemented for names in both SpatVector and SpatRaster
- `rast,data.frame` method to avoid confusion with the `matrix` and `list` methods in response to a [SO question](https://stackoverflow.com/q/68133958/635245) by Stackbeans
- the extreme values used to represent NA where not as intended (one or two lower) for INT2U and INT4U. Reported by Jean-Luc Dupouey on [stackoverflow](https://stackoverflow.com/q/68216362/635245)
- `writeCDF` now also writes the time dimensions if there is only one time-step. See this [SO question](https://stackoverflow.com/a/68227180/635245)
- `vect,character` (filename) now has argument `layer` to select a layer from a multi-layer file / database, and arguments `query`, `extent` and `filter` for reading a subset
- `subst` can now create multiple output layers See [issue 276](https://github.com/rspatial/terra/issues/276) by Agustin Lobo
- `classify` can now create different multiple output layers See [issue 276](https://github.com/rspatial/terra/issues/276) by Agustin Lobo
- Argument `alpha` of `plot,SpatRaster` can now be a `SpatRaster`. See this [SO question](https://stackoverflow.com/q/68736432/635245) by James McCarthy


## bug fixes 

- The `filename` and `overwrite` arguments were ignored in `rasterize`
- gdal options are now also honored for create-copy drivers [#260](https://github.com/rspatial/terra/issues/260)
- buffer for lonlat now works better at the world's "edges" [#261](https://github.com/rspatial/terra/issues/261)
- scale/offset were ignored by `project`. Reported by Fabian Fischer
- `rasterize(SpatRaster,SpatVector)` with `inverse=TRUE` crashed the R session. Issue [#264](https://github.com/rspatial/terra/issues/264) by Jean-Luc Dupouey
- The output of `merge` and `mosaic` was not correct for large rasters (only the first rows were used). Reported by Zavud Baghirov in [#271](https://github.com/rspatial/terra/issues/271)
- `as.points,SpatRaster` did not remove `NA`'s correctly and shifted values. Issues [#269](https://github.com/rspatial/terra/issues/269) and [#273](https://github.com/rspatial/terra/issues/273) by Julian Hagenauer
- `rast-matrix` rotated values when using an equal-sided matrix [#274](https://github.com/rspatial/terra/issues/274) by Jakub Nowosad
- the number of rows and columns were reversed when using `project` with a crs argument. [#283](https://github.com/rspatial/terra/issues/283) by Timothée Giraud
- In `classify`, argument `right` had TRUE and FALSE reversed. 
- `terrain` had edge effects [#303](https://github.com/rspatial/terra/issues/303) by Andrew Gene Brown.
- `terrain` can now compute multiple variables at once [#286](https://github.com/rspatial/terra/issues/286) by Žan Kuralt
- `wrap,SpatRaster` changed factors into numeric [#302](https://github.com/rspatial/terra/issues/302) by Patrick Schratz
- `writeVector` failed with "FlatGeobuf" (and probably other formats as well) for not using a proper MultiPolygon [#299](https://github.com/rspatial/terra/issues/299) by L Dalby
- regular sampling of polygons with `spatSample` is now much more regular [#289](https://github.com/rspatial/terra/issues/289) by Jakub Nowosad



# version 1.3-4

## new

- `na.omit,SpatVector-method` to remove empty geometries and/or attribute records that have an `NA`
- new method `src` to create a `SpatRasterCollection` (a loose collection of tiles). 
- `merge` and `mosaic` now have methods for a `SpatRasterCollection`. To avoid the (inefficient) use of `do.call`. See issue [#210](https://github.com/rspatial/terra/issues/210) by Matthew Talluto.
- `activeCat` and `activeCat<-` to get or set the "active" category if there are multiple categories (raster attributes)
- `as.numeric` and `catalyze` to transfer categories to numeric cell values
- summarize methods such as `range` and `mean` for (the attributes of) a `SpatVector`
- new method `shade`, to compute hill shading

## enhancements

- additional arguments (such as `na.rm`) are now used by `rasterize` with point geometries. Suggested by Jakub Nowosad  [#209](https://github.com/rspatial/terra/issues/209)
- improved handling (and documentation) of `gstat` models by `interpolate`. See issue [#208](https://github.com/rspatial/terra/issues/208) by Jakub Nowosad.
- new argument `cpkgs` to `predict` to list the packages that need to be exported to the cores if argument `cores` is larger than one. `?predict` now shows different approaches to parallelize `predict` (based on examples in issue [#178](
https://github.com/rspatial/terra/issues/178) raised by by Matthew Coghill).
- `freq` now returns labels for categorical layers
- `adjacent` now has a `pairs` argument. Requested by Kenneth Blake Vernon in issue [#239](https://github.com/rspatial/terra/issues/239) 
- `adjacent` now also takes a matrix to specify adjacent cells
- `mean` and other summarize methods now take a `filename` argument and disallow non-recognized named arguments. See issue [#238](https://github.com/rspatial/terra/issues/238) by Jessica Nephin
- The raster attribute table of ESRI-GRID integer data, or from an ESRI `vat.dbf` file is now ignored if it only has the counts of the values. See issue [#234]( https://github.com/rspatial/terra/issues/234) by Jullee
- time attributes are no longer lost when doing raster operations. Suggested by Mauricio Zambrano-Bigiarini in [#246]( https://github.com/rspatial/terra/issues/246)
- resample (and project) no longer ignore `gdal=""` write options and use BIGTIFF if necessary (suggested by Ani Ghosh)
- new argument `layer` in the `extract-SpatRaster,SpatVector` method to extract values for a single layers specified for each geometry (see this [question](https://gis.stackexchange.com/a/401591/8993)).

## bug fixes 

- better handling of paths with non-ASCII characters (e.g., Chinese) for GeoTiff but still fails for NetCDF (see issue [#233](https://github.com/rspatial/terra/issues/223) by Dongdong Kong)
- `extract` with points and `cells=TRUE` or `xy=TRUE` gave garbled output
- `as.character,SpatRaster-method` (called by `wrap`) did not capture the layer names. Reported by Pascal Title [#213](https://github.com/rspatial/terra/issues/213)
- `focal` mirrored the weight matrix, thus affecting the results when using an asymmetrical weight matrix. Reported by Sebastiano Trevisani
- `terra::terraOptions` now works without attaching the package (issue [#229](https://github.com/rspatial/terra/issues/229) reported by Karl Dunkle Werner)
- `app` with `ncores > 0` and a function that returns multiple layers now works (issue [#240](https://github.com/rspatial/terra/issues/240) reported by BastienFR.
- `autocor` (local) can now handle `NA` values. Reported by Jakub Nowosad [#245](https://github.com/rspatial/terra/issues/245).
- `mask` with a SpatVector and a large (out of memory) multi-layer SpatRaster only worked for the first layer. Reported by Monika Tomaszewska.



# version 1.2-10

## new

- `as.lines` method for SpatRaster
- `as.polygons` method for SpatVector lines
- `autocor,numeric-method` has new methods `mean`, to compute the local mean, and `locmor`, for the local Moran's *I* 
- `sharedPaths` method for SpatVector (lines and polygons)
- `RGB2col` method to reduce a three-layer RGB SpatRaster to a single layer SpatRaster with a color-table (with <= 256 colors)
- `split` methods for SpatVector and SpatRaster


## enhancements

- `rast(Raster*)` now takes the crs from the Raster object, not from the file it may point to. Suggested by Floris Vanderhaeghe [#200](https://github.com/rspatial/terra/issues/200)
- `convhull` has a new argument `by=""` to make convex hulls for sub-sets of a SpatVector.
- faster processing of large in memory rasters. See issue [#206](https://github.com/rspatial/terra/issues/206) by Krzysztof Dyba.


## bug fixes

- `extract` with multiple layers could return a data.frame where the values were not in the correct order (by row instead of by column)
- `crop` works again with `sf` objects. Reported by Sebastian Brinkmann [#201](https://github.com/rspatial/terra/issues/201)
- `vect,sf-method` now also works for lines, and should be faster
- `vect,character` crashed R if a file had empty geometries. Reported by consumere [#202](https://github.com/rspatial/terra/issues/202)
- `extract(points, bilinear=TRUE, cells=TRUE)` now works. Reported by fab4app [#203](https://github.com/rspatial/terra/issues/203)
- `zonal` now works for `min` and `max`. Reported by Jakub Nowosad  [#207](https://github.com/rspatial/terra/issues/207)


## name changes

To avoid name conflicts with the `spatstat` package

- `area,SpatRaster-method(x, sum=FALSE)` -> `cellSize(x)`
- `area,SpatRaster/SpatVector-method(x, sum=TRUE)` -> `expanse(x)`
- `convexhull` -> `convHull`
- `perimeter` -> `perim`
- `tiles` -> `makeTiles`
- `coords` -> `crds`


# version 1.2-5

## new

- `trim` has a new argument `value` that allows trimming rows and columns with other values than the default `NA`
- `rapp` has a new argument `clamp` that allows clamping start and end values to `1:nlyr(x)`, avoiding that all values are considered `NA`
- `spatSample,SpatRaster-method` has new arguments `as.points` and `values`. Getting values, cells and coordinates is no longer mutually exclusive. In response to [#191](https://github.com/rspatial/terra/issues/191). Requested by Agustin Lobo
- `area,SpatRaster-method` has a new argument `mask=FALSE`
- `classify` can now take a single number to request that many cuts
- `mosaic` and `merge` now warn and resample if rasters are not aligned
- `extract` has a new argument `exact` to get the fraction covered for each cell

## bug fixes

- `flip(x, direction="vertical")` no longer reverses the order of the layers
- `extract` did not work for horizontal or vertical lines as their extent was considered invalid. Reported by Monika Tomaszewska
- `autocor` did not handle NA values [#192](https://github.com/rspatial/terra/issues/192). Reported by Laurence Hawker
- `nearest` now works for angular coordinates
- The unit of `slope` in `terrain` was not correct (the tangent was returned instead of the slope) [#196](https://github.com/rspatial/terra/issues/196). Reported by Sven Alder
- `quantile` now works for rasters that have cells that are all `NA`. Reported by Jerry Nelson

## name changes

To avoid name conflicts with tidyverse 

with deprecation warning:

- separate -> segregate
- expand -> extend
- near -> nearby
- pack -> wrap 

without deprecation warning:

- transpose -> trans
- collapse -> tighten 
- fill -> fillHoles
- select -> sel


# version 1.1-17

## major changes 

- `c,SpatVector-method` now returns a list. `rbind` is used to append SpatVector objects
- overhaul of handling of factors. `rats` has been removed, and `levels` and `cats` have changed


# version 1.1-4

- No news recorded for this version and earlier versions

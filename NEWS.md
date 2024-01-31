# version 1.7-71

## bug fixes

- k_means did not work if there were NAs [#1314](https://github.com/rspatial/terra/issues/1314) by Jakub Nowosad
- `layerCor` with a custom function did not work anymore [#1387](https://github.com/rspatial/terra/issues/1387) by Jakub Nowosad
- `plet` broke when using "panel=TRUE" [#1384](https://github.com/rspatial/terra/issues/1384) by Elise Hellwig
- using /vis3/ to open a SpatRaster did not work [#1382](https://github.com/rspatial/terra/issues/1382) by Mike Koontz
- `plot<SpatRaster>(add=TRUE)` sampled the raster data without considering the extent of the map. [#1394](https://github.com/rspatial/terra/issues/1394) by Márcia Barbosa
- `plot<SpatRaster>(add=TRUE)` now only considers the first layer of a multi-layer SpatRaster [1395](https://github.com/rspatial/terra/issues/1395) by Márcia Barbosa
- `set.cats` failed with a tibble was used instead of a data.frame [#1406](https://github.com/rspatial/terra/issues/1406) by Mike Koontz
- `polys` argument "alpha" was ignored if a single color was used. [#1413](https://github.com/rspatial/terra/issues/1413) by Derek Friend
- `query` ignore the "vars" argument if all rows were selected. [#1398](https://github.com/rspatial/terra/issues/1398) by erkent-carb.
- `spatSample` ignored "replace=TRUE" with random sampling, na.rm=TRUE, and a sample size larger than the non NA cells. [#1411](https://github.com/rspatial/terra/issues/1411) by Babak Naimi
- `spatSample` sometimes returned fewer values than requested and available for lonlat rasters. [#1396](https://github.com/rspatial/terra/issues/1396) by Márcia Barbosa.


## enhancements

- `vect<character>` now has argument "opts" for GDAL open options, e.g. to declare a file encoding. [#1389](https://github.com/rspatial/terra/issues/1389) by Mats Blomqvist
- `plot(plg=list(tic=""))` now allows chosing alternative continuous legend tic-mark styles ("in", "out", "through" or "none")
- `makeTiles` has new argument "buffer" [#1408](https://github.com/rspatial/terra/issues/1408) by Joy Flowers.


## new 

- `prcomp<SpatRaster>` method [#1361](https://github.com/rspatial/terra/issues/1361#issuecomment-1860311029) by Jakub Nowosad
- `add_box` to add a box around the map. The box is drawn where the axes are, not around the plotting region.
- `getTileExtents` provides the extents of for tiles. These may be used in parallelization. See [#1391](https://github.com/rspatial/terra/issues/1391) by Alex Ilich.


# version 1.7-65

Released 2023-12-15

## bug fixes

- `flip` with argument `direction="vertical"` filed in some cases with large rasters processed in chunks [0b714b0](https://github.com/rspatial/terra/commit/0b714b038b101011ce27391133a7a1f9a91821cc) by Dulci on [stackoveflow](
https://stackoverflow.com/questions/77304534/rspatial-terraflip-error-when-flipping-a-multi-layer-spatrast-object) 
- SpatRaster now correctly handles `NA & FALSE` and `NA | TRUE` [#1316](https://github.com/rspatial/terra/issues/1316) by John Baums
- `set.names` wasn't working properly for SpatRasterDataset or SpatRasterCollection [#1333](https://github.com/rspatial/terra/pull/1333) by Derek Friend
- `extract` with argument "layer" not NULL shifted the layers [#1332](https://github.com/rspatial/terra/issues/1332) by Ewan Wakefield
- `terraOptions` did not capture "memmin" on [stackoverflow](https://stackoverflow.com/questions/77552234/controlling-chunk-size-in-terra) by dww
- `rasterize` with points and a built-in function could crash if no field was used [#1369](https://github.com/rspatial/terra/issues/1369) by anjelinejeline


## enhancements

- `mosaic` can now use `fun="modal"`
- `rast<matrix> and rast<data.frame>` now have option 'type="xylz" [#1318](https://github.com/rspatial/terra/issues/1318) by Agustin Lobo
- `extract<SpatRaster,SpatVector>` can now use multiple summarizing functions [#1335](https://github.com/rspatial/terra/issues/1335) by Derek Friend
- `disagg` and `focal` have more optimistic memory requirement estimation [#1334](https://github.com/rspatial/terra/issues/1334) by Mikko Kuronen

## new

- `k_means<SpatRaster>` method [#1314](https://github.com/rspatial/terra/issues/1314) by Agustin Lobo
- `princomp<SpatRaster>` method [#1361](https://github.com/rspatial/terra/issues/1361) by Alex Ilich
- `has.time<SpatRaster>` method 
- new argument "raw=FALSE" to `rast`, `sds`, and `sprc` to allow ignoring scale and offset [1354](https://github.com/rspatial/terra/issues/1354) by Insang Song


# version 1.7-55

Released 2023-10-14

## bug fixes

- `mosaic` ignored the filename argument if the SpatRasterCollection only had a single SpatRaster [#1267](https://github.com/rspatial/terra/issues/1267) by Michael Mahoney
- Attempting to use `extract` with a raster file that had been deleted crashed R. [#1268](https://github.com/rspatial/terra/issues/1268) by Derek Friend
- `split<SpatVector,SpatVector>` did not work well in all cases. [#1256](https://github.com/rspatial/terra/issues/1256) by Derek Corcoran Barrios
- `intersect` with two SpatVectors crashed R if there was a date/time variable [#1273](
https://github.com/rspatial/terra/issues/1273) by Dave Dixon
- "values=FALSE" was ignored by `spatSample<SpatRaster>(method="weights")` [#1275](https://github.com/rspatial/terra/issues/1275) by François Rousseu
- `coltab<-` again works with a list as value [#1280](https://github.com/rspatial/terra/issues/1280) by Diego
Hernangómez
- `stretch` with histogram equalization was not memory-safe [#1305](https://github.com/rspatial/terra/issues/1305) by Evan Hersh
- `plot` now resets the "mar" parameter [#1297](https://github.com/rspatial/terra/issues/1297) by Márcia Barbosa
- `plotRGB` ignored the "smooth" argument [#1307](https://github.com/rspatial/terra/issues/1307) by Timothée Giraud


## enhancements

- argument "gdal" in `project` was renamed to "use_gdal" [#1269](https://github.com/rspatial/terra/issues/1269) by Stuart Brown.
- SpatVector attributes can now be stored as an ordered factor [#1277](https://github.com/rspatial/terra/issues/1277) by Ben Notkin
- `plot<SpatVector>` now uses an "interval" legend when breaks are supplied [#1303](https://github.com/rspatial/terra/issues/1303) by Gonzalo Rizzo
- `crop<SpatRaster>` now keeps more metadata, including variable names [#1302](https://github.com/rspatial/terra/issues/1302) by rhgof
- `extract(fun="table")` now returns an easier to use data.frame [#1294] (https://github.com/rspatial/terra/issues/1294) by Fernando Aramburu.


## new
- `metags<-` and `metags` to set arbitrary SpatRaster/file level metadata [#1304](https://github.com/rspatial/terra/issues/1304) by Francesco Chianucci 


# version 1.7-46

Released 2023-09-06

## bug fixes

- `plot<SpatVector>` used the wrong main label in some cases [#1210](https://github.com/rspatial/terra/issues/1210) by Márcia Barbosa
- `plotRGB` failed with an "ext=" argument [#1228](https://github.com/rspatial/terra/issues/1228) by Dave Edge
- `rast<array>` failed badly when the array had less than three dimensions. [#1254](https://github.com/rspatial/terra/issues/1254) by andreimirt.
- `all.equal` for a SpatRaster with multiple layers [#1236](https://github.com/rspatial/terra/issues/1236) by Sarah Endicott
- `zonal(wide=FALSE)` could give wrong results if the zonal SpatRaster had "layer" as layername. [#1251](https://github.com/rspatial/terra/issues/1251) by Jeff Hanson
- `panel` now support argument "range" [#141](https://github.com/rspatial/terra/issues/1241) by Jakub Nowosad
- `rasterize` with `by=` returned wrong layernames if the by field was not sorted [#1266](https://github.com/rspatial/terra/issues/1266) by Sebastian Dunnett
- `mosaic` with multiple layers was not correct [#1262](https://github.com/rspatial/terra/issues/1262) by Jean-Romain


## enhancements

- `wrap<SpatRaster>` now stores color tables [#1215](https://github.com/rspatial/terra/issues/1215) by Patrick Brown
- `global` now has a "maxcell" argument [#1213](https://github.com/rspatial/terra/issues/1213) by Alex Ilich
- `layerCor` with fun='pearson' now returns output with the layer names [#1206](https://github.com/rspatial/terra/issues/1206)
- `vrt` now has argument "set_names" [#1244](https://github.com/rspatial/terra/issues/1244) by sam-a-levy
- `vrt` now has argument "return_filename" [#1258](https://github.com/rspatial/terra/issues/1258) by Krzysztof Dyba
- `project<SpatRaster>` has new argument "by_util" exposing the GDAL warp utility [#1222](https://github.com/rspatial/terra/pull/1222) by Michael Sumner.


## new
- `compareGeom` for list and SpatRasterCollection [#1207](https://github.com/rspatial/terra/issues/1207) by Sarah Endicott
- `is.rotated<SpatRaster>` method [#1229](https://github.com/rspatial/terra/issues/1229) by Andy Lyons
- `forceCCW<SpatVector>` method to force counter-clockwise orientation of polygons [#1249](https://github.com/rspatial/terra/issues/1249) by srfall.
- `vrt_tiles` returns the filenames of the tiles in a vrt file [#1261](https://github.com/rspatial/terra/issues/1261) by Derek Friend
- `extractAlong` to extract raster cell values for a line that are ordered along the line. [#1257](https://github.com/rspatial/terra/issues/1257) by adamkc.


# version 1.7-39

Released 2023-06-23

## bug fixes

- the tempdir option did not use path.expand. [#1195](https://github.com/rspatial/terra/issues/1195) by Alex Ilich
- the layer names returned by predict where inconsistent when using argument "index". [#1194](https://github.com/rspatial/terra/issues/1194) by Michael Mahoney
- compilation failed with older compilers because of use of std::filesystem [#1191](https://github.com/rspatial/terra/issues/1191)
- Small changes to `RGB<-` and `coltab<-` so that terra can be installed with R-devel (after a bug fix https://bugs.r-project.org/show_bug.cgi?id=18538)


# version 1.7-37

Released 2023-06-18

## bug fixes

- `rasterize` with points and a custom function did not work for large rasters. [#1127](https://github.com/rspatial/terra/issues/1127) by Skip Woolley
- `crop<SpatRaster, SpatVector>` with "mask=TRUE" did not work well if the raster had a scale/offset [#1128](https://github.com/rspatial/terra/issues/1128) by Monika Anna Tomaszewska
- `zonal<SpatRaster>` with a custom function always removed NAs. [#1133](https://github.com/rspatial/terra/issues/1133) by  Matthias Weigand
- `wrap<SpatRaster>` lost changed layer names if the source was from disk; and information on some time-step in some cases. [#1144](https://github.com/rspatial/terra/issues/1144) by Pascal Führlich
- `global(fun="isNA")` was not correct when the SpatRaster had multiple layers [#1141](https://github.com/rspatial/terra/issues/1141) by Robin Freeman
- `interpIDW` with `near=TRUE` did not work properly (near=TRUE is now the default). [#1186](https://github.com/rspatial/terra/issues/1186) by Hugh Graham
- "YYYY-1-1" was sometimes encoded as "YYYY-13-1". [#1168](https://github.com/rspatial/terra/issues/1168) by Colin Brust


## enhancements

- `panel` for categorical SpatRasters. [#1143](https://github.com/rspatial/terra/issues/1143) by Jason Flower
- argument "ext" in `plot<SpatRaster>` can now also expand the plot. [#1136](https://github.com/rspatial/terra/issues/1136) by Jakub Nowosad.
- argument `overwrite=FALSE` to `makeTiles`. [#1167](https://github.com/rspatial/terra/issues/1167) by Gray Martin.
- legend options for `<plet,SpatVector`>. [#1177](https://github.com/rspatial/terra/issues/1177) by Agustin Lobo.
- better handling of mixed geometry type vector data by `vect` and `svc`. [#1160](https://github.com/rspatial/terra/issues/1160) by Mike Sumner.
- new argument `sql` to `query<SpatVectorProxy>`. [#1157](https://github.com/rspatial/terra/issues/1157) by Carl Boettiger
- support for writing raster data to a vitual file system [#1209](https://github.com/rspatial/terra/issues/1209) by Carl Boettiger

## new
- `wrap<SpatRasterDataset>` and `wrap<SpatRasterCollection>` methods. [#954](https://github.com/rspatial/terra/issues/954) by James Camac


# version 1.7-29

Released 2023-04-22

## new

- `regress<SpatRaster,numeric>` to get regression model coefficients for each cell, with a fixed "X".
- `regress<SpatRaster,SpatRaster>` to get regression model coefficients for each cell.

## enhancements

- `lapp<SpatRasterDataset>` is now more flexible in that it can now also use functions that are vectorized by cell, not by chunk. See [#1029](https://github.com/rspatial/terra/issues/1029)
- `project<SpatVector>` has new argument "partial=FALSE" that can be used to keep geometries that can only be partially included in the output crs.
- extracting a SpatVector column with a non-existing variable name now returns NULL (because that is what a data.frame does) instead of throwing an error. [#1118](https://github.com/rspatial/terra/issues/1118) by Derek Friend.

## bug fixes

- a problem with reading empty categories in .img files created buggy SpatRasters
- `global` with fun="notNA" was wrong [#111](https://github.com/rspatial/terra/issues/1111) by Jeffrey Hanson
- `extract<SpatRaster,SpatVector>` with "bind=TRUE" did not work
- `extract<SpatRaster,SpatVector>` with point geometries and a "fun" returned values in the wrong order
- `plot<SpatRaster>` argument "colNA" did not work when "alpha" was also set [#1102](https://github.com/rspatial/terra/issues/1102) by Márcia Barbosa
- `crop<SpatRaster>` with "extend=TRUE" did not extend the SpatRaster if the input had no cell values. [#1114](https://github.com/rspatial/terra/issues/1114) by Jasper van Doninck
- setting a factor or date/time variable in a SpatVector did not work [#1117](https://github.com/rspatial/terra/issues/1117) by MK Schneider
- `focalMat` did not work well when using terraOptions(todisk=T) [#1116](https://github.com/rspatial/terra/issues/1116)


# version 1.7-23

Released 2023-04-08

## new

- The `halo` function for adding halo-ed text to plots is now exposed
- `add_legend` to allow using a keyword such as "topleft" to position a custom legend. [#1053](https://github.com/rspatial/terra/issues/1053) by Márcia Barbosa
- the `same.crs` function is now exported
- `countNA<SpatRaster>` method
- `split<SpatVector,SpatVector>` to split polygons with lines

## enhancements

- better support for other color spaces than RGB [#1060](https://github.com/rspatial/terra/issues/1060) by Dominic Royé
- path expansion in writeVector [#1055](https://github.com/rspatial/terra/issues/1055) by Andrew Gene Brown.
- `clamp<SpatRaster>` now also accepts SpatRasters to set the lower and upper boundaries.
- `freq` has new arguments "zones=NULL" and "wide=FALSE", to allow tabulation of values by zone.
- `expanse<SpatRaster>` has new arguments "zones=NULL" and "wide=FALSE", to allow tabulation of values by zone.
- `unique<SpatRaster>` has new argument "digits=NA"
- `rasterize<SpatRaster,SpatVector>` now accepts fun="table" to tabulate cells by cell value
- `rast<character>` has new argument "snap" to snap the window in or out. [#1094](https://github.com/rspatial/terra/issues/1094) by Derek Friend
- `plot` has new argument "clip=TRUE" that can be set to FALSE to avoid clipping the axes to the mapped area [#1080](https://github.com/rspatial/terra/issues/1080) by Márcia Barbosa
- better error message when coercing an sf object that is not fully formed [#1098](https://github.com/rspatial/terra/issues/1098) by Brandon McNellis
- `writeCDF<SpatRaster>` had new argument "split" allowing to treat each layer as a subdataset [#1077](https://github.com/rspatial/terra/issues/1077) by Andrea Manica
- `global` now accepts multiple summarizing functions 

## bug fixes

- A SpatRaster with RGB layers was forced to INT1U when writing [#1051](https://github.com/rspatial/terra/issues/1051) by Cesar Aybar
- In files with multiple vector layers, the crs of the first layer was always used; ignoring that the crs could be different for other layers [#1052](https://github.com/rspatial/terra/issues/1052) by Andrew Gene Brown
- `sieve` was not able to write to file [#1061](https://github.com/rspatial/terra/issues/1061) by leo
- `rasterize` did not work with sf objects [#1054](https://github.com/rspatial/terra/issues/1054) by Jakub Nowosad
- `query` did not work for hyphenated layer names [#1058](https://github.com/rspatial/terra/issues/1058) by Robbie Price
- `focal3D` na.policy did not work [#1057](https://github.com/rspatial/terra/issues/1057) by Flávio Mota
- `layerCor` with `na.rm=TRUE` failed for a SpatRaster with more than 2 layers [#1056](https://github.com/rspatial/terra/issues/1056) by Alex Ilich.
- inset with keyword positioning did not work well [#1053](https://github.com/rspatial/terra/issues/1053) by Márcia Barbosa
- yearmonths time stamps were not read from file for years <1970 and >2037 [#1062](https://github.com/rspatial/terra/issues/1062) by Colin Brust
- `compareGeom` did not work for multiple SpatRasters [#1063](https://github.com/rspatial/terra/issues/1064)
- `viewshed` could not handle a filename argument. [#1100](https://github.com/rspatial/terra/issues/1100) by kamistick


# version 1.7-18

Released 2023-03-06

## new

- argument `order=FALSE` to `sort<SpatRaster>` 
- `sort<SpatVector>` (and `<data.frame>` method
- argument `by=NULL` to `rasterize>` [#986](https://github.com/rspatial/terra/issues/986) by Sam Weber
- `meta<SpatRaster>` method to get metadata
- `compare<SpatRaster>` and `logic<SpatRaster>` methods
- `vect<SpatExtent>` method
- `panel<SpatRaster>` for "panel" plots (multiple layers, single legend)

## enhancements

- it is now possible to save terra options across sessions [#995](https://github.com/rspatial/terra/issues/995) by Guillaume Patoine.
- better warnings for `is.lonlat` [#1006](https://github.com/rspatial/terra/issues/1006) by Andrew Gene Brown
- argument `na.rm` to `merge<SpatRaster>`
- the axes of maps created with `plot` are now snug around the mapped area, instead of at the limits of the graphics figure region.
- C++ cleaning to avoid warnings by clang-tidy (e.g. now using `.empty()` instead of `.size()==0`). [#1013-1017] by Michael Chirico 
- `rasterize` with lines and polygons can now use the "fun" argument (for min, max, mean, and sum) [#1041](https://github.com/rspatial/terra/issues/1041) by Bart Huntley

## bug fixes

- the legend created by `plet` was not always correct. [#983](https://github.com/rspatial/terra/issues/983) by Simon Rolph
- `spatSample<SpatRaster>(regular=TRUE)` failed with providing two numbers (row, col) as sample size. [#991](
https://github.com/rspatial/terra/issues/991) by srfall
- `merge<SpatRaster>` did not ignore NAs [#1002](https://github.com/rspatial/terra/issues/1002) by jmmonnet.
- `writeCDF` failed when using argument force_v4 [#1009](https://github.com/rspatial/terra/issues/1009) by R. Kyle Bocinsky
- `predict` better handling of rasters with many NAs [#988](https://github.com/rspatial/terra/issues/998) by Lucas Johnson
- `layerCor` did not handle NAs well if they were in different cells across layers [#1034](https://github.com/rspatial/terra/issues/1034) by François Rousseu.


# version 1.7-3

Released 2023-01-24

## new

- argument `w` to `zonal<SpatRaster,SpatRaster>` to compute weighted means
- `zonal<SpatRaster,SpatVector>` method
- `clamp_ts` method

## bug fixes 

- in the previous version, a bug was introduced such that the order of operation in arithmetic operations with SpatRasters was ignored. [#978](https://github.com/rspatial/terra/issues/978) by Andrew Marx
- Fixed `split<SpatVector>`. [#979](https://github.com/rspatial/terra/issues/979) by srfall
- `spatSample` with `as.df=FALSE` returned a data.frame instead of a matrix [#982](https://github.com/rspatial/terra/issues/982) by Alex Ilich


# version 1.6-53

Released 2023-01-17

## new

- arithmetic and logical operations between a SpatRaster and a matrix, to allow for using cell-varying and cell/layer-varying scalars. layer-varying scalars were already supported via vectors.

## enhancements

- `shade` is now vectorized for arguments `angle` and `direction` to facilitate generating multiple hillshades that can be combined for a better result [#948](https://github.com/rspatial/terra/issues/948) by Jürgen Niedballa
- `sharedPaths` now uses spatial indices [#960](https://github.com/rspatial/terra/issues/960) by Jeff Hanson
- `predict` has better support for models such as ranger that do not return anything for missing values [#968](https://github.com/rspatial/terra/issues/968) by Alex Ilich

## bug fixes 

- `writeCDF` now supports writing yearly time steps [#926](https://github.com/rspatial/terra/issues/926) by Andrea Manica
- `as.contour` now works for a single level [#966](https://github.com/rspatial/terra/issues/966) by Johannes Signer
- subsetting a SpatRaster with a window returned a SpatRaster with the dimensions of the non-windowed raster, thus changing the resolution. [#964](https://github.com/rspatial/terra/issues/964) by Derek Friend
- removing a factor variable from a SpatVector crashed R. [#969](https://github.com/rspatial/terra/issues/969) by Andrew Gene Brown
- median did not always return the correct number for a SpatRaster with 3 or more layers [#970](https://github.com/rspatial/terra/issues/970) by MatteaE


# version 1.6-47

Released 2022-12-02

## new

- `roll<SpatRaster>` method for rolling (moving) average and other rolling functions
- `noNA<SpatRaster>` method to identify cells that are not NA (across layers) 
- `rangeFill<SpatRaster>` method 

## enhancements

- argument `exhaustive` to `spatSample<SpatRaster>` for large sparse rasters. [#905] by PetiteTong.
- `focalPairs` and `focalReg` can now use the values in custom windows as weights. [#907] by Fabian Fischer.
- `focalReg` now has additional argument "intercept=TRUE". [#916] by Jordan Adamson
- `crs(x, warn=TRUE)<-` now emits a warning about the difference between transforming and setting a crs when x already had a crs. [#897] by Márcia Barbosa.
- it is now possible to write a scale and offset with `writeRaster` [#900] by Kyle David
- `crosstab` now shows the labels names for a categorical SpatRaster. [895] by Derek Corcoran Barrios
- `makeTiles` can now take a SpatVector to define the tiles. [920] by Tristan Goodbody

## bug fixes 

- `focalPairs` and `focalReg` now work for custom windows [#907] by Fabian Fischer
- argument "alpha" in `plot<SpatVector>` was not working properly. [#906] by Márcia Barbosa.
- `time<-` with time-step "years" could not handle negative years. [#911] by Andrea Manica
- `wrap`/`unwrap` (and by extension `saveRDS`/`readRDS`) did not handle categorical rasters well [#912] by Christine Anderson.
- `interpIDW` failed with GDAL 3.6 [#910] by Roger Bivand
- `spatSample` with strata bug fix "unable to find an inherited method for function 'trim'" [#919] by Alfredo Ascanio
- it is possible to slice a SpatRaster with a SpatExtent [#914] by Jakub Nowosad.
- `merge`/`mosaic` did not handle NAs when using two layers [#913] by Joao Carreiras.


## name changes

- focalCor -> focalPairs to reflect its possible use beyond correlation


# version 1.6-41

Released 2022-11-18

## new

- `[` and `[<-` for SpatRaster now have a third index `k` for subsetting or assigning values by layer
- `anyNA` and `allNA` for SpatRaster
- `unwrap` to restore a PackedSpatVector or PackedSpatRaster
- `rasterizeWin` method for rasterization with a moving window (circle, ellipse, rectangle, buffer)
- `interpIDW` method for inverse-distance-weighted interpolation of points with a moving window
- `interpNear` method for nearest neighbor interpolation of points with a moving window
- `viewshed` method for SpatRaster
- `update` method for SpatRaster to write new names or a new extent or crs to an existing raster file.
- `sieve` filter for SpatRaster
- argument `segments=FALSE` to `disagg<SpatRaster>`
- `sprc<character>` method to create a SpatRasterCollection from a file with subdatasets
- `graticule` function to create a SpatGraticule and related methods `plot<SpatGraticule>` and `crop<SpatGraticule>`
- `elongate` method for SpatVector lines


## enhancements

- faster `mosaic` and `merge<SpatRaster>` [#577] by Jean-Romain
- `wrap<SpatRaster>` now uses file references if the data is deemed to large to all load into memory. [#801] by Jean-Romain
- `readRDS` and `unserialize` now return a SpatRaster or SpatVector (instead of a PackedSpat*)
- better support for a "local" arbitrary Euclidean crs [#797] by Agustin Lobo
- `clamp` can now take low and high values for each layer 
- The `pax` argument in `plot` now provides more control over what to draw on each axis via parameters `side`, `tick` and `lab`
- The `pax` argument in `plot` now has argument `retro` to use a sexagesimal notation of degrees
- `extend` has a new argument `fill=NA`
- A warning is now given when `c`ombining SpatRasters with different CRSs. [#818] by Andrew Marx
- `plotRGB` now accounts for the value of zlim when stretching; allowing to use the same coloring scheme across SpatRasters [#810] by Agustin Lobo.
- the center of rotation for `spin` is now vectorized


## bug fixes 

- The annoying garbage collection messages `Error in x$.self$finalize() : attempt to apply non-function` is now suppressed in most cases. [#218] by Charlie Joey Hadley. This problem should go away altogether when a new version of "Rcpp" is released (ETA Jan 2023) thanks to a fix by Kevin Ushey [#30]
- `spatSample` with `na.rm` and SpatRasters with multiple layers did not work. [#800] by Andrea Manica
- `adjacent<SpatRaster>` with `pairs=TRUE, include=TRUE` ignored `include=TRUE` [#808] by Joseph Lewis
- `rasterize` did not accept "NA" as value for updating [#809]  by Márcia Barbosa
- `extract` with a perfectly vertical or horizontal line failed in some cases [#823] by Dimitri Falk
- `wrap<SpatVector>` failed if there was a single point geometry [#815] by Patrick Schaefer
- `extract<SpatRaster>` with `weights=TRUE` did not return values [#814] by Jean-Luc Dupouey. 
- `x[["newname"]] <- r` for SpatRasters `x` and `r` did not work [#795] by Jim Shady
- fixed support for some non-conventional netCDF files [#869] by Mike Sumner, [#864] by eleanorecc, and [#851] by Philippe Massicotte.

## name changes

- `costDistance` -> `costDist` to avoid conflict with {gdistance}
- `gridDistance` -> `gridDist` for consistency


# version 1.6-17

Released 2022-09-10

## new

- `droplevels` for SpatRaster. [#757] by Rodolfo Jaffe.
- `normalize.longitude` for SpatVector. 
- `scoff` to get and `scoff<-` to set the scale (gain) and offset of a SpatRaster. 

## enhancements

- new argument `raw=FALSE` to `extract<SpatRaster>` [#776] by Thomas Roh.
- `as.data.frame` now takes `na.rm=NA` to only remove rows that are NA for all layers. The default value changed from `TRUE` to `NA`. [#792] by Ed Carnell
- faster plotting of SpatVector data [#774] by Krzysztof Dyba
- `distance<SpatRaster>` has new arguments "target" and "exclude". [#560] by Bernardo Brandão Niebuhr
- new argument `sparse=FALSE` for `relate<SpatVector,SpatVector>. 
- new argument `usenames=FALSE` for `lapp<SpatRasterDataset>` [#793] by Colin Brust.
- `vect<character>` now reports that a file is non-existent [#784] by John Baums
- faster `relate` [#716] by Krzysztof Dyba
- `focal3D` now checks if all the window's dimensions are odd [#772] by Neander Marcel Heming

## bug fixes 

- `all.equal` bug [#756] fixed by John Baums
- `extract<"SpatRaster","sf">` ignored the ID argument. [#755] by Dainius Masiliūnas.
- There is now (in all cases) a check to avoid overwriting (one of) the input file(s) when writing a raster file [#760] by John Baums
- `vrt` is no longer constrained by the maximum number of files that can be opened [#780] by 8Ginette8	
- `weighted.mean` crashed with numeric weights and na.rm=TRUE [#777] by David Holstius
- `project<SpatRaster>` did not consider an extent that was set by the user [#775] by Philippe Massicotte
- `focalCor` failed for large rasters [#607] by John Clark
- `focal` with `expand=TRUE` was prone to run out of memory [#610] by Nathan Elliott
- `crop<SpatVector>` did not work well when the second argument were points or lines [#782] by Márcia Barbosa
- `adjacent` with `pairs=TRUE` now respects the `include=TRUE` argument [808] by Joseph Lewis


# version 1.6-7

Released 2022-08-07

## new

- method `blocks` to guide reading raster data in chunks. [#748] by John Baums

## enhancements 

- A warning is given when writing raster values that are outside the limits of the requested datatype [#752] by Jim Shady
- Arguments to `extract` were simplified. [#736] by François Rousseu

## bug fixes 

- values of `focal` where not correct if the input SpatRaster had multiple layers and a "custom" function. [#727] by Jean-Luc Dupouey. 
- `plot<SpatRaster>` did not honor argument `legend=FALSE`. [#738] by Grzegorz Sapijaszko
- `expanse` failed when processing in chunks [#741] by Gareth Davies 
- `crop<SpatRaster,SpatExtent>` with argument `snap="out"` could lead to a crash if the extent was beyond the SpatRaster. [#740] by Mauricio Zambrano-Bigiarini


# version 1.6-3

Released 2022-07-25

## bug fixes

- `subst` no longer uses values that it changed earlier on. [#639] by Paul Smith
- `as.points<SpatRaster>` could return wrong factor labels. [#640] by Attilio Benini
- `mask<SpatRaster,SpatVector>` crashed when the results were written to disk. [#646] by Monika Anna Tomaszewska
- `extract<SpatRaster,SpatVector(points)>(xy=TRUE)` returned the locations of the points, not the xy-coordinates of the cells. [#650] by Ward Fonteyn
- `wrap<SpatRaster>` did not return the correct labels for some categorical rasters. [#652] by Jakub Nowosad
- better support for non-latin characters in the legend [#658] by Krzysztof Dyba
- holes in small lon/lat polygons are now properly buffered [#689] by David Hofmann

## enhancements 

- `subst` can now substitute the values from multiple input layers with a single output value (layer)
- `subset<SpatVector>` now behaves like `subset<data.frame>` [#648] by Andrew Gene Brown
- setting category labels with a vector of names is now deprecated. A data.frame with at least two columns should be used. The first column should have the cell values (IDs).
- It is now possible to "drop" a layer from a SpatRaster by setting it to NULL [#664] by Daniel Valentins
- `freq` now provides the labels of factors, even if `bylayer=FALSE`. It now always returns a `data.frame` (it used to return a `matrix` in some cases. [#687] by Rodolfo Jaffé
- `disagg` and `aggregate` now return a warning instead of an error when using a (dis)aggregation factor of 1.[#684] by Justin Fain.
- `project` crashed when erroneously projecting raster data from one celestial body to another [#688] by Mike Sumner
- you can now set a color table with a two column (value, ID) data.frame
- categorical rasters can now be updated more easily [#667] by Alex Ilich
- more control over matching values with colors when using `plot`. [#673] by Jakub Nowosad.
- SpatVector attributes can now also be a factor, date, or POSIXct. [#697] by Grant Williamson
- improved handling of missing values in `extract(method="bilinear")`. [#693] by swooping-magpie

## new

- argument `as.raster` to `unique<SpatRaster>` to create a categorical raster with the unique combinations in the layers of the input raster. The default for argument `na.rm` was changed to `FALSE`
- `sort<SpatRaster>` to sort cell values across layers.
- `has.colors` and `has.RGB` for SpatRaster
- `cover` can now combine categorical rasters 
- `concats` to combine the levels of two SpatRaster into new categories [#663] by Alex Ilich
- `zonal<SpatVector,SpatVector>` method to aggregate SpatVector attributes by polygons


# version 1.5-34

Released 2022-06-09

## bug fixes

- "flipped" rasters were not always handled well. [#546] by Dan Baston 
- better reading of GTiff with subdatsets. [#601] by Kyle Doherty
- better handling of multi-layer categorical rasters and `extract`. [#580] by André M. Bellvé
- `weighted.mean` did not adjust the weights if there were `NA`s in the values. [#574] by Lars Dalby
- bug in masking. [#552] reported by Márcia Barbosa and [565] by Jakub Nowosad.
- fixed `stretch` option in `plotRGB` [#550] by Agustin Lobo
- unwrap of a SpatRaster failed with a crs including a "'". [#602] by Jean Romain.
- `spatSample` with `cells=TRUE` failed for planar data [#544] by Benjamin Misiuk
- `compareGeom(x, y, stopOnError=FALSE)` did not remove the error messages stored in `x` leading to unexpected warnings later on. [#568] by David Hofmann.

## enhancements 

- Using & or | with SpatRasters now returns a boolean SpatRaster. [#594] by Dan Baston 
- SpatVector now supports logical values. [#593] by Derek Friend
- Attempt to create SpatRaster with an invalid number of rows now gives an error. [#544] by Dan Baston
- `layerCor` does not create temp files anymore. [#551] by Christine Anderson
- not using the same iterator symbols in nested loops to avoid warnings on the Intel compiler. [#573] by Gareth Davies.

## new

- new arguments `res` and `origin` to `project<SpatRaster>` method. [#596] by Alex Ilich
- new argument `inside=TRUE` to `centroids` to get a centroid-like point that is guaranteed to be on the geometry ("point on surface"). [#588] by Márcia Barbosa
- new argument `keepgeom=FALSE` to `vect<data.frame>` that allows setting (keeping) the geometry as an attribute. [#586] by Márcia Barbosa
- `saveRDS` and `serialize` methods for SpatRaster and SpatVector. [#549] by Andrei Mîrț
- `xFromCol` and `yFromCol` now have a `<SpatRaster,missing>` method. [#583] by Michael Sumner.
- `svc<sf>` method to deal with GeometryCollection types. [#585] by Sarah Endicott
- `as.points<SpatRaster>` and `as.polygons<SpatRaster>` have a new argument `na.all=FALSE` that affects the interpretation of `na.rm`. [#548] by Jean-Luc Dupouey.
- `setGDALconfig` and `getGDALconfig` to set GDAL configuration options. [#608] by Erik Bolch.
- new argument `circular` to `rapp` to allow the start to be after the end (for if layers represent days of the year)
- new method `costDistance<SpatRaster>` 
- new methods `where.min` and `where.max` for `SpatRaster` to get the cell numbers for the extreme values in a SpatRaster. 
- new method `emptyGeoms<SpatVector>` to get the indices of empty (null) geometries
- new method `rasterizeGeom` to rasterize the geometry count or the area of (small) polygons or the length of lines.
- new method `not.na` for `SpatRaster` which is a shortcut for `!is.na(x)`.
- `as.list` implemented for `<SpatRasterDataset>`.
- `sources` implemented for `<SpatRasterDataset>`, `<SpatVector>` and `<SpatVectorProxy>` [#638] by Andrew Gene Brown

## name changes

- delauny -> delaunay [#627] by Derek Friend


# version 1.5-21

Released 2022-02-17

- `writeVector` and `vect` now work with GPGK if the path has non-ascii characters [#518]
- The results of `predict` with `cores > 1` and more than one output variable were garbled
- `zonal` dropped category names when using an external (R) function [#527] by Jakub Nowosad
- focal/focalCpp showed strange patterns when the window size was larger than the block size [#519] by Alex Ilich
- using `xy=TRUE` in `as.data.frame` normalized the names [#538] by Kodi Arfer
- new argument `options` to `vrt` [#629] by Monika Tomaszewska.

## enhancements 

- `makeTiles` has new arguments `extend` and `na.rm` [#520] by by L. Dalby
- `project<SpatRaster>` now uses nearest neighbor as default method for RGB rasters
- new argument `na.rm=TRUE` to `unique`. [#561] by Matthieu Stigler


# version 1.5-17

Released 2022-01-30

## bug fixes

- `app<SpatRasterDataset>` ignored the filename. [#498] by jszhao
- `vect<data.frame>` failed silently if xy coordinates were integers [#496] by Márcia Barbosa
- The output of `aggregate<SpatRaster>` was malformed when `nrow(x) %% fact != 0`. [#492] by Jean-François Bourdon
- Integer `NA`s in SpatVector attributes where only recognized on Windows [#491] by Márcia Barbosa
- `plot<SpatVector>` failed when using a character variable with many unique values. [#489] by Márcia Barbosa
- `rotate` failed on large files. Reported by Ujjawal Singh
- writing raster files with a color table could lead to a crash [#501] by Kodi Arfer
- `crds` replicated the coordinates [#504] by Murray Efford
- `as.data.frame<SpatRaster>` returned integers if the file stored values as integers, even if there was a scale/offset that creates decimal numbers [#509] by Kodi Arfer
- `project` opened the input raster file in read/write mode intead of read mode. That did not work with files that cannot be updated.

## enhancements 

- `distance`, `gridDistance`, `direction` and `patches` now process all layers of the input SpatRaster. [#503] by Chris Haak
- consistent copy-on-modify behavior in `()<-` methods. in-place updating available with `set.` methods such as `set.names` and `set.values`. [#493] by Jean Romain and [#511] by Bryan Fuentes
- much faster writing of GPGK vector data by using a single transaction (following sf) [#460] by Krzysztof Dyba
- `aggregate<SpatRaster>` now accepts functions that return more than one value per aggregated cell
- `writeVector` has new argument `insert` to add a layer to an existing file (e.g. GPKG).

## new

- new option `method="weights"` for `spatSample<SpatRaster>`
- new `mask<SpatVector,SpatVector>` method to select intersecting geometries
- new method `is.related`
- `values<SpatRaster>` has new option `na.rm=TRUE`. [#490] by Henk Harmsen
- new class `SpatVectorProxy` to provide access to large vector databases that cannot or should not be read into memory in its entirety.
- new argument `proxy=FALSE` to `vect` to create a SpatVectorProxy object
- new method `query<SpatVectorProxy>` to extract parts of a SpatVectorProxy
- new method `vector_layers` that returns, and can delete, vector format layers from a database/file such as GPKG

## name changes

To avoid name clashes with tidyverse 

- arrow -> north
- src -> sprc
- simplify -> simplifyGeom

For consistency 

- setCats -> set.cats


# version 1.5-12

Released 2022-01-13

## bug fixes

- `setValues` and `init` failed (or even crashed R) when using a single value on a largish raster. [#414]
- conversion from `sfc` to `SpatVector` lost the crs. [#415] by Jean-Luc Dupouey
- `buffer` on a SpatRaster with no values caused a crash [#416] by Sebastian Brinkmann
- `writeVector` now assumes "traditional GIS order" (long/lat) if the CRS specifies lat/long. [#333](
 by Agustin Lobo
- argument `main` was ignored in `density` when using a single layer SpatRaster [#424] by dvictori
- Summary type math functions such as `min` and `mean`, when used with multiple SpatRasters and numbers, ignored additional SpatRasters [#426] by Zhuonan Wang
- names are now conserved when creating a SpatRaster from a RasterStack that points to file(s) [#430] by Dan Baston
- `classify` with `right=FALSE` ignored `include.lowest=TRUE` [#442] by Alex Ilich
- `patches` now combines patches that connect across the data line [#366] by Hirscht
- `patches(directions=8)` now connects in NE/SW direction [#451] by Jean-François Bourdon.
- `centroids` now considers cases where SpatVector parts are nearest to each other when crossing the date line in stead of the zero-meridian [#366] by Hirscht 
- `terrain` created empty (`NA`) rows between chunks used for processing large rasters. [#453] by Robert Ritson.
- `inset` did not draw the "box" correctly. [#457] by Márcia Barbosa
- `as.lines` now works with a points SpatVector [#464] by Márcia Barbosa 


## enhancements 

- `values(x)<-` now accepts (hex coded) colors as values
- `focal` now wraps around the dateline like raster::focal [#242] by Alexander Marbler
- `aggregate` now does not show a progress bar in all cases [#249] by Lachlan
- `as.data.frame<SpatRaster> or <SpatVector>` are now also implemented as S3 methods to assure correct dispatch by other S3 methods such as `data.table::as.data.table`. See [#284] by Patrick Schratz
- `crs` now shows the correct authority if it is not EPSG. [#419] by Matthew Williamson
- It now possible to add a SpatRaster to an empty SpatRaster (with no values), even if it has a different geometry, ignoring the empty SpatRaster [#421] by Alex Ilich.
- `rast<filename>` has a new argument `lyrs` to subset the layers and open the file in one step.
- `rast<array>` now has a crs and extent argument. [#439] by RS-eco
- `type="xyz"` is now default in `rast<data.frame>`. [#438] by RS-eco
- `classify` has a new argument `brackets` to show if a side of an interval is open or closed.
- further support for categorical data in `freq` and `as.data.frame`. [#441] ngould7
- speed up in processing of multi-layer in memory data. [#437] by Krzysztof Dyba
- `vect<matrix>` and `vect<data.frame>` are now much faster. [#413] by BastienFR 	
- `extract` with points provided as a matrix or cell numbers is not much faster. [#341]
- `focal` has a new argument `na.policy` that can be set to one of "all" (default), "only" or "omit". argument `na.only` has been removed, as you can now use `na.policy="only"`
- `inset` argument `border` changed to `perimeter` to allow passing `border` on to `plot<Spat*>`. [#456] by Márcia Barbosa
- The compile-time and run-time versions of GEOS are now compared and a warning is given if they are not the same. [#459] by Edzer Pebesma
- it is now possible to add sub-datasets to GPKG and GTiff files. [#300] by gtitov
- general option `memfrac` can now be set to zero (in stead of not lower than 0.1). [#476] by Matt Strimas-Mackey
- new argument `allowGaps` in `patches` to disallow gaps between patch IDs. See [#478] by Dunbar Carpenter.


## new 

- timestamps and units are now saved to an auxiliary file (filename.aux.json) for all raster formats except NetCDF when using writeCDF (because in that case they are stored in the netcdf file)
- new method `mergeTime` to combine multiple rasters, perhaps partly overlapping in time, into a single time series
- new method `fillTime` that can add empty layers in between existing layers to assure that the time step between layers is constant 
- new method `approximate` to fill in missing values by cell across layers
- new methods `is.bool` and `as.bool` for SpatRaster and explicit recognition of Boolean raster data in various places (e.g., extract, plot)
- new methods `is.int` and `as.int` for SpatRaster. 
- when assigning integer values to a SpatRaster, or when reading an integer file, the corresponding layers are now classified as being of integer type [#446] by L. Dalby
- new method `layerCor` (like `raster::layerStats`). [#420] by Alex Ilich
- new method `focalCor` (like `raster::corLocal`). [#427] by Zhuonan Wang
- new method `all.equal` for `SpatRaster`. See [#428] by Dongdong Kong
- new method `math` for `SpatRaster` that implements the Math-generic methods *and* accepts a filename
- new method `sds<array>` 
- new method `rasterize<matrix>`, see [#413] by BastienFR 	
- new method `colorize` to transform color representations 	
- new method `arrow` to draw a (North) arrow on a map. [#461] by Márcia Barbosa
- new method `densify` to insert nodes between existing nodes of a line or polygon SpatVector
- new method `direction` for SpatRaster. [#462] by Márcia Barbosa
- new method `focal3D` to compute focal values for a three-dimensional (row, column, layer) window
- new function `makeVRT` to create a vrt file for a file that needs a header to be read.
- new option `method="stratified"` for `spatSample<SpatRaster>`. [#470] by Michael Mahoney
- new general option `memmax` to cap the amount of RAM that terra can be used in raster processing [#476] by Matt Strimas-Mackey
- new method `gridDistance` to compute distances traversing a raster, perhaps with obstacles. [#477] by Márcia Barbosa


# version 1.4-22

Released 2021-11-24

## changes 
- `focal` now has ellipses (`...`) to allow for providing additional arguments to `fun`. For this reason it does not have a `na.rm` argument anymore as that can be supplied via the ellipses. In practice this means that the default will be `na.rm=FALSE` for the standard functions such as `mean` and `sum`.


## bug fixes

- `app` grossly overestimated RAM needed, slowing it down. Reported by Jerry Nelson 
- `terra` now installs, again, with older versions of GEOS [#406] by fparyani
- `terra` did not install with Clang on CRAN/OSX due to using C++13 idiom.


## enhancements 

- `lapp` and `tapp` now have a `cores` argument (as do `app` and `predict`). Suggested by Dongdong Kong [#365]
- `focal` now also works with a function that returns multiple values. See [#318] by Alex Ilich. 
- `focal` can now process multiple layers in one step. 
- expanded support for conversion from `stars` objects [#220] by Jakub Nowosad


## new 

- `focalCpp` takes a C++ function that iterates over cells to speed up computation by avoiding `apply` (see [#318] by Alex Ilich). 
- `focalReg` for focal OLS regression models between layers 


# version 1.4-20

Released 2021-11-16

## bug fixes

- `terra` did not install with versions of GDAL below 3 [#402] by Alex Ilich.
- `distance` between two SpatVectors or matrices with `pairwise=FALSE` returned a matrix that was filled by column instead of by row [#403] by Paul Smith


# version 1.4-19

Released 2021-11-15

## bug fixes

- `rast` with some NetCDF files failed because of bad parsing of dates. [#361] by Juan Carlos Zamora-Pereira
- `distance<SpatRaster>` with lon/lat data was not correct. [#368]
by Greg Schmidt
- `as.polygons<SpatRaster>` failed with a SpatRaster and a categorical layer that is not the first layer. [#370] by Patrick Schratz
- The filename argument in `rasterize` was not ignored, also causing errors when writing to temporary files. [#377] by Robbie Price
- `rast<character>` crashed if the sds was an empty character string. [#381] by Dan Baston
- `plot<SpatVector>` now responds to the `range` argument [#385] by Márcia Barbosa
- `zonal` failed for user-defined functions. [#393] by mqueinnec


## new

- new method `selectHighest` to select n cell values with the highest or lowest values. 
- new method `vect<list>` to append SpatVectors (faster than `do.call(rbind, x)`)
- new argument `align=FALSE` to `project` to align to the template SpatRaster but ignore the resolution
- new method `gdalCache` to set the GDAL cache size, contributed by Dan Baston [#387]
- new method `fileBlocksize`
- new argument `options` to `writeVector` to pass layer creation options to GDAL
- new SpatVector topology methods `mergeLines`, `snap`, `makeNodes`, `removeDupNodes`, `gaps`, `simplify`
- new SpatVector characterization methods `width` and `clearance`


## enhancements 

- `terra` now installs with older versions of GEOS [#363]
- `terra` now installs on CentOS 7 with GDAL 2.1.4 and a C++ compiler that does not support std::regexp. [#384] by Ariel Paulson


# version 1.4-11

Released 2021-10-11

## enhancements

- the definition of `setValues` now has two arguments (`x` and `values`), just like `raster` had; to avoid reverse dependency problems with `raster`

# version 1.4-9

Released 2021-10-07

## name changes

To avoid name conflicts with `sp` (via `raster`) `disaggregate` is now called `disagg` and `bbox,SpatRaster` and `bbox<SpatVector>` have been removed (but could be resurrected in `raster` or under another name).

## enhancements

- `project` and `resample` now choose the resampling method based on the first layer, using "near" for categorical data. Thanks to Matthew Lewis [#355]

## bug fixes

- `hist` failed with small samples. Issue [#356] by Martin Queinnec


# version 1.4-7

Released 2021-10-05

## note

`terra` no longer depends on `raster`. To avoid name clashes between these two packages, and to allow replacing methods from `rgeos` and `rgdal` in `raster`, `raster` now depends on `terra` instead. 


## enhancements

- `freq` has a new argument `usenames`. See issue [#309] by Bappa Das
- `rast<character>` has a new argument `opts` that can be used to pass GDAL open options. See issue [#314]
- `rast<SpatRaster>` now takes arguments `names` and `vals`. See issue [#323] by Dongdong Kong
- `crs<-` now warns if an unsupported datum is used. See issue [#317]
- `spatSample` now returns factor values if a SpatRaster layer is.factor except when using `as.df=FALSE`
- new method `origin<-` to set the origin of a SpatRaster. See issue [#326] by Jakub Nowosad
- `crs` has a new argument `parse`. See [#344]
- `plot<SpatRaster,missing>` has a new argument `reset=FALSE` that allows resetting the par()$mar parameters after plotting. See issue [#340] by Derek Friend
- `crds` has a new argument `na.rm`. See [#338] by Kodi Arfer 
- `show(Spat*)` now prints the name and EPSG code of a crs if available. See [#317] by Jakub Nowosad


## bug fixes 

- `plotRGB` failed if there were `NA`s. Issue [#308] by Jakub Nowosad
- `writeVector` crashed R when used with a SpatVector with no geometries. Reported by Timothy White in issue [#319]
- `summary<SpatRaster>` now returns counts for the classes (instead of a numerical summary of the indices) [#324] by Jakub Nowosad
- `tapp` with a character index now returns a SpatRaster with the correct names [#345] by Stuart Brown 
- `rasterize` with a character variable now adds the ID column to the categories [#337] by Tate Brasel
- `cellSize` now masks values in all cases (when requested with `mask=TRUE`). Issue [#339] by Jean-Luc Dupouey
- `buffer<SpatVector>` no longer treats lines like polygons [#332] by Márcia Barbosa
- `plot` now passes the layer index to `fun` [#310] by Ben Tupper
- the `to_id` in `nearest` was sometimes wrong. See [#328] by Shawn Ligocki
- better support for ESRI value attribute tables (VAT). See this [SO question]( https://stackoverflow.com/q/69385928/635245)
- `focal` did not reset initial values for NA cells when processing chunks. [#312] by Jeffrey Evans
- `focal` could run out of memory when using a large window and user-defined function, and was inexact at the chunk boundary [#347]
- `zonal` with `as.raster=TRUE` failed for categorical SpatRasters [#348] by Jakub Nowosad



# version 1.3-22

Released 2021-08-20

## enhancements

- if `time(x) <- d` is set with a `Date` class object, `time(x)` now returns a `Date` object instead of a `POSIXct` object. Issue [#256] by Mauricio Zambrano-Bigiarini
- The UTF-8 encoding of character attributes of a SpatVector is now declared such that they display correctly in R. See issue [#258] by AGeographer. Also implemented for names in both SpatVector and SpatRaster
- `rast<data.frame>` method to avoid confusion with the `matrix` and `list` methods in response to a [SO question](https://stackoverflow.com/q/68133958/635245) by Stackbeans
- the extreme values used to represent NA where not as intended (one or two lower) for INT2U and INT4U. Reported by Jean-Luc Dupouey on [stackoverflow](https://stackoverflow.com/q/68216362/635245)
- `writeCDF` now also writes the time dimensions if there is only one time-step. See this [SO question](https://stackoverflow.com/a/68227180/635245)
- `vect<character>` (filename) now has argument `layer` to select a layer from a multi-layer file / database, and arguments `query`, `extent` and `filter` for reading a subset
- `subst` can now create multiple output layers See [issue 276] by Agustin Lobo
- `classify` can now create different multiple output layers See [issue 276] by Agustin Lobo
- Argument `alpha` of `plot<SpatRaster>` can now be a `SpatRaster`. See this [SO question](https://stackoverflow.com/q/68736432/635245) by James McCarthy


## bug fixes 

- The `filename` and `overwrite` arguments were ignored in `rasterize`
- gdal options are now also honored for create-copy drivers [#260]
- buffer for lonlat now works better at the world's "edges" [#261]
- scale/offset were ignored by `project`. Reported by Fabian Fischer
- `rasterize<SpatRaster,SpatVector>` with `inverse=TRUE` crashed the R session. Issue [#264] by Jean-Luc Dupouey
- The output of `merge` and `mosaic` was not correct for large rasters (only the first rows were used). Reported by Zavud Baghirov in [#271]
- `as.points,SpatRaster` did not remove `NA`'s correctly and shifted values. Issues [#269] and [#273] by Julian Hagenauer
- `rast<matrix>` rotated values when using an equal-sided matrix [#274] by Jakub Nowosad
- the number of rows and columns were reversed when using `project` with a crs argument. [#283] by Timothée Giraud
- In `classify`, argument `right` had TRUE and FALSE reversed. 
- `terrain` had edge effects [#303] by Andrew Gene Brown.
- `terrain` can now compute multiple variables at once [#286] by Žan Kuralt
- `wrap<SpatRaster>` changed factors into numeric [#302] by Patrick Schratz
- `writeVector` failed with "FlatGeobuf" (and probably other formats as well) for not using a proper MultiPolygon [#299] by L Dalby
- regular sampling of polygons with `spatSample` is now much more regular [#289] by Jakub Nowosad



# version 1.3-4

Released 2021-06-20

## new

- `na.omit<SpatVector>` to remove empty geometries and/or attribute records that have an `NA`
- new method `src` to create a `SpatRasterCollection` (a loose collection of tiles). 
- `merge` and `mosaic` now have methods for a `SpatRasterCollection`. To avoid the (inefficient) use of `do.call`. #210 by Matthew Talluto.
- `activeCat` and `activeCat<-` to get or set the "active" category if there are multiple categories (raster attributes)
- `as.numeric` and `catalyze` to transfer categories to numeric cell values
- summarize methods such as `range` and `mean` for (the attributes of) a `SpatVector`
- new method `shade`, to compute hill shading

## enhancements

- additional arguments (such as `na.rm`) are now used by `rasterize` with point geometries. #209 by Jakub Nowosad
- improved handling (and documentation) of `gstat` models by `interpolate`. #208 by Jakub Nowosad
- new argument `cpkgs` to `predict` to list the packages that need to be exported to the cores if argument `cores` is larger than one. `?predict` now shows different approaches to parallelize `predict` (based on examples in issue. #178 by by Matthew Coghill.
- `freq` now returns labels for categorical layers
- `adjacent` now has a `pairs` argument. #239 by Kenneth Blake Vernon
- `adjacent` now also takes a matrix to specify adjacent cells
- `mean` and other summarize methods now take a `filename` argument and disallow non-recognized named arguments. #238 by Jessica Nephin
- The raster attribute table of ESRI-GRID integer data, or from an ESRI `vat.dbf` file is now ignored if it only has the counts of the values. #234 by Jullee
- time attributes are no longer lost when doing raster operations. #246 by Mauricio Zambrano-Bigiarini
- resample (and project) no longer ignore `gdal=""` write options and use BIGTIFF if necessary (suggested by Ani Ghosh)
- new argument `layer` in the `extract-SpatRaster,SpatVector` method to extract values for a single layers specified for each geometry (see this [question](https://gis.stackexchange.com/a/401591/8993)).

## bug fixes 

- better handling of paths with non-ASCII characters (e.g., Chinese) for GeoTiff but still fails for NetCDF. [#233] by Dongdong Kong
- `extract` with points and `cells=TRUE` or `xy=TRUE` gave garbled output
- `as.character<SpatRaster>` (called by `wrap`) did not capture the layer names. [#213] by Pascal Title
- `focal` mirrored the weight matrix, thus affecting the results when using an asymmetrical weight matrix. Reported by Sebastiano Trevisani
- `terra::terraOptions` now works without attaching the package. [#229] by Karl Dunkle Werner
- `app` with `ncores > 0` and a function that returns multiple layers now works. [#240] by BastienFR.
- `autocor` (local) can now handle `NA` values. [#245] by Jakub Nowosad .
- `mask` with a SpatVector and a large (out of memory) multi-layer SpatRaster only worked for the first layer. Reported by Monika Tomaszewska.



# version 1.2-10

Released 2021-05-13

## new

- `as.lines` method for SpatRaster
- `as.polygons` method for SpatVector lines
- `autocor<numeric>` has new methods `mean`, to compute the local mean, and `locmor`, for the local Moran's *I* 
- `sharedPaths` method for SpatVector (lines and polygons)
- `RGB2col` method to reduce a three-layer RGB SpatRaster to a single layer SpatRaster with a color-table (with <= 256 colors)
- `split` methods for SpatVector and SpatRaster


## enhancements

- `rast<Raster*>` now takes the crs from the Raster object, not from the file it may point to. [#200] by Floris Vanderhaeghe 
- `convhull` has a new argument `by=""` to make convex hulls for sub-sets of a SpatVector.
- faster processing of large in memory rasters. See issue [#206] by Krzysztof Dyba.


## bug fixes

- `extract` with multiple layers could return a data.frame where the values were not in the correct order (by row instead of by column)
- `crop` works again with `sf` objects. [#201] by Sebastian Brinkmann 
- `vect<sf>` now also works for lines, and should be faster
- `vect<character>` crashed R if a file had empty geometries. [#202] by consumere
- `extract(points, bilinear=TRUE, cells=TRUE)` now works. [#203] by fab4app 
- `zonal` now works for `min` and `max`. [#207] Reported by Jakub Nowosad


## name changes

To avoid name conflicts with the `spatstat` package

- `area,SpatRaster-method(x, sum=FALSE)` -> `cellSize(x)`
- `area,SpatRaster/SpatVector-method(x, sum=TRUE)` -> `expanse(x)`
- `convexhull` -> `convHull`
- `perimeter` -> `perim`
- `tiles` -> `makeTiles`
- `coords` -> `crds`


# version 1.2-5

Released 2021-04-30

## new

- `trim` has a new argument `value` that allows trimming rows and columns with other values than the default `NA`
- `rapp` has a new argument `clamp` that allows clamping start and end values to `1:nlyr(x)`, avoiding that all values are considered `NA`
- `spatSample<SpatRaster>` has new arguments `as.points` and `values`. Getting values, cells and coordinates is no longer mutually exclusive. In response to [#191] by Agustin Lobo
- `area<SpatRaster>` has a new argument `mask=FALSE`
- `classify` can now take a single number to request that many cuts
- `mosaic` and `merge` now warn and resample if rasters are not aligned
- `extract` has a new argument `exact` to get the fraction covered for each cell

## bug fixes

- `flip(x, direction="vertical")` no longer reverses the order of the layers
- `extract` did not work for horizontal or vertical lines as their extent was considered invalid. Reported by Monika Tomaszewska
- `autocor` did not handle NA values [#192] by Laurence Hawker
- `nearest` now works for angular coordinates
- The unit of `slope` in `terrain` was not correct (the tangent was returned instead of the slope), [#196] by Sven Alder
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

Released 2021-04-14

## major changes 

- `c<SpatVector>` now returns a list. `rbind` is used to append SpatVector objects
- overhaul of handling of factors. `rats` has been removed, and `levels` and `cats` have changed


# version 1.1-4

- No news recorded for this and earlier versions

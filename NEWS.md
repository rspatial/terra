# version 1.6-15

## new

- `droplevels` for SpatRaster. [#757](https://github.com/rspatial/terra/issues/757) by Rodolfo Jaffe.
- `normalize.longitude` for SpatVector. 
- `scoff` to get and `scoff<-` to set the scale (gain) and offset of a SpatRaster. 

## enhancements

- new argument `raw=FALSE` to `extract<SpatRaster>` [#776](https://github.com/rspatial/terra/issues/776) by Thomas Roh.
- `as.data.frame` now takes `na.rm=NA` to only remove rows that are NA for all layers. The default value changed from `TRUE` to `NA`. [#792](https://github.com/rspatial/terra/issues/792) by Ed Carnell
- faster plotting of SpatVector data [#774](https://github.com/rspatial/terra/issues/774) by Krzysztof Dyba
- `distance<SpatRaster>` has new arguments "target" and "exclude". [#560](https://github.com/rspatial/terra/issues/560) by Bernardo Brandão Niebuhr
- new argument `sparse=FALSE` for `relate<SpatVector,SpatVector>. 
- new argument `usenames=FALSE` for `lapp<SpatRasterDataset>` [#793](https://github.com/rspatial/terra/issues/793) by Colin Brust.
- `vect<character>` now reports that a file is non-existent [#784](https://github.com/rspatial/terra/issues/784) by John Baums
- faster `relate` [#716](https://github.com/rspatial/terra/issues/716) by Krzysztof Dyba
- `focal3D` now checks if all the window's dimensions are odd [#772](https://github.com/rspatial/terra/issues/772) by Neander Marcel Heming

## bug fixes 

- all.equal bug [#756](https://github.com/rspatial/terra/pull/756) fixed by John Baums
- extract<"SpatRaster","sf"> ignored the ID argument. [#755](https://github.com/rspatial/terra/issues/755) by Dainius Masiliūnas.
- There is now (in all cases) a check to avoid overwriting (one of) the input file(s) when writing a raster file [#760](https://github.com/rspatial/terra/issues/760) by John Baums
- `vrt` is no longer constrained by the maximum number of files that can be opened [#780](https://github.com/rspatial/terra/issues/780) by 8Ginette8	
- `weighted.mean` crashed with numeric weights and na.rm=TRUE [#777](https://github.com/rspatial/terra/issues/777) by David Holstius
- `project<SpatRaster>` did not consider an extent that was set by the user [#775](https://github.com/rspatial/terra/issues/775) by Philippe Massicotte
- `focalCor` failed for large rasters [#607](https://github.com/rspatial/terra/issues/607) by John Clark
- `focal` with `expand=TRUE` was prone to run out of memory [#610](https://github.com/rspatial/terra/issues/610) by Nathan Elliott
- `crop<SpatVector>` did not work well when the second argument were points or lines [#782](https://github.com/rspatial/terra/issues/782) by Márcia Barbosa


# version 1.6-7

Released on 2022-08-07

## new

- method `blocks` to guide reading raster data in chunks. [#748](https://github.com/rspatial/terra/issues/748) by John Baums

## enhancements 

- A warning is given when writing raster values that are outside the limits of the requested datatype [#752](
https://github.com/rspatial/terra/issues/752) by Jim Shady
- Arguments to `extract` were simplified. [#736](https://github.com/rspatial/terra/issues/736) by François Rousseu

## bug fixes 

- values of `focal` where not correct if the input SpatRaster had multiple layers and a "custom" function. [#727](https://github.com/rspatial/terra/issues/727) by Jean-Luc Dupouey. 
- `plot<SpatRaster>` did not honor argument `legend=FALSE`. [#738](https://github.com/rspatial/terra/issues/738) by Grzegorz Sapijaszko
- `expanse` failed when processing in chunks [#741](https://github.com/rspatial/terra/issues/741) by Gareth Davies 
- `crop<SpatRaster,SpatExtent>` with argument `snap="out"` could lead to a crash if the extent was beyond the SpatRaster. [#740](https://github.com/rspatial/terra/issues/740) by Mauricio Zambrano-Bigiarini



# version 1.6-3

Released on 2022-07-25

## bug fixes

- `subst` no longer uses values that it changed earlier on. [#639](https://github.com/rspatial/terra/issues/639) by Paul Smith
- `as.points<SpatRaster>` could return wrong factor labels. [#640](https://github.com/rspatial/terra/issues/640) by Attilio Benini
- `mask<SpatRaster,SpatVector>` crashed when the results were written to disk. [#646](https://github.com/rspatial/terra/issues/646) by Monika Anna Tomaszewska
- `extract<SpatRaster,SpatVector(points)>(xy=TRUE)` returned the locations of the points, not the xy-coordinates of the cells. [#650](https://github.com/rspatial/terra/issues/650) by Ward Fonteyn
- `wrap<SpatRaster>` did not return the correct labels for some categorical rasters. [#652](https://github.com/rspatial/terra/issues/652) by Jakub Nowosad
- better support for non-latin characters in the legend [#658](https://github.com/rspatial/terra/issues/658) by Krzysztof Dyba
- holes in small lon/lat polygons are now properly buffered [#689](https://github.com/rspatial/terra/issues/689) by David Hofmann


## enhancements 

- `subst` can now substitute the values from multiple input layers with a single output value (layer)
- `subset<SpatVector>` now behaves like `subset<data.frame>` [#648](https://github.com/rspatial/terra/issues/648) by Andrew Gene Brown
- setting category labels with a vector of names is now deprecated. A data.frame with at least two columns should be used. The first column should have the cell values (IDs).
- It is now possible to "drop" a layer from a SpatRaster by setting it to NULL [#664](
https://github.com/rspatial/terra/issues/664) by Daniel Valentins
- `freq` now provides the labels of factors, even if `bylayer=FALSE`. It now always returns a `data.frame` (it used to return a `matrix` in some cases. [#687](https://github.com/rspatial/terra/issues/687) by Rodolfo Jaffé
- `disagg` and `aggregate` now return a warning instead of an error when using a (dis)aggregation factor of 1.[#684](https://github.com/rspatial/terra/issues/684) by Justin Fain.
- `project` crashed when erroneously projecting raster data from one celestial body to another [#688](https://github.com/rspatial/terra/issues/688) by Mike Sumner
- you can now set a color table with a two column (value, ID) data.frame
- categorical rasters can now be updated more easily [#667](https://github.com/rspatial/terra/issues/667) by Alex Ilich
- more control over matching values with colors when using `plot`. [673](https://github.com/rspatial/terra/issues/673) by Jakub Nowosad.
- SpatVector attributes can now also be a factor, date, or POSIXct. [697](https://github.com/rspatial/terra/issues/697) by Grant Williamson
- improved handling of missing values in `extract(method="bilinear")`. [693](https://github.com/rspatial/terra/pull/693) by swooping-magpie

## new

- argument `as.raster` to `unique<SpatRaster>` to create a categorical raster with the unique combinations in the layers of the input raster. The default for argument `na.rm` was changed to `FALSE`
- `sort<SpatRaster>` to sort cell values across layers.
- `has.colors` and `has.RGB` for SpatRaster
- `cover` can now combine categorical rasters 
- `concats` to combine the levels of two SpatRaster into new categories [663](https://github.com/rspatial/terra/issues/663) by Alex Ilich
- `zonal<SpatVector,SpatVector>` method to aggregate SpatVector attributes by polygons


# version 1.5-34

Released on 2022-06-09

## bug fixes

- "flipped" rasters were not always handled well. [#546](https://github.com/rspatial/terra/issues/546) by Dan Baston 
- better reading of GTiff with subdatsets. [#601](https://github.com/rspatial/terra/issues/601) by Kyle Doherty
- better handling of multi-layer categorical rasters and `extract`. [#580](https://github.com/rspatial/terra/issues/580) by André M. Bellvé
- `weighted.mean` did not adjust the weights if there were `NA`s in the values. [#574](https://github.com/rspatial/terra/issues/574) by Lars Dalby
- bug in masking. [#552](https://github.com/rspatial/terra/issues/552) reported by Márcia Barbosa and [565](https://github.com/rspatial/terra/issues/565) by Jakub Nowosad.
- fixed `stretch` option in `plotRGB` [#550](https://github.com/rspatial/terra/issues/550) by Agustin Lobo
- unwrap of a SpatRaster failed with a crs including a "'". [#602](https://github.com/rspatial/terra/issues/602) by Jean Romain.
- `spatSample` with `cells=TRUE` failed for planar data [#544](https://github.com/rspatial/terra/issues/554) by Benjamin Misiuk
- `compareGeom(x, y, stopOnError=FALSE)` did not remove the error messages stored in `x` leading to unexpected warnings later on. [#568](https://github.com/rspatial/terra/issues/568) by David Hofmann.

## enhancements 

- Using & or | with SpatRasters now returns a boolean SpatRaster. [#594](https://github.com/rspatial/terra/issues/594) by Dan Baston 
- SpatVector now supports logical values. [#593](https://github.com/rspatial/terra/issues/593) by Derek Friend
- Attempt to create SpatRaster with an invalid number of rows now gives an error. [#544](https://github.com/rspatial/terra/issues/544) by Dan Baston
- `layerCor` does not create temp files anymore. [#551](https://github.com/rspatial/terra/issues/551) by Christine Anderson
- not using the same iterator symbols in nested loops to avoid warnings on the Intel compiler. [#573](https://github.com/rspatial/terra/issues/573) by Gareth Davies.

## new

- new arguments `res` and `origin` to `project<SpatRaster>` method. [#596](https://github.com/rspatial/terra/issues/596) by Alex Ilich
- new argument `inside=TRUE` to `centroids` to get a centroid-like point that is guaranteed to be on the geometry ("point on surface"). [#588](https://github.com/rspatial/terra/issues/588) by Márcia Barbosa
- new argument `keepgeom=FALSE` to `vect<data.frame>` that allows setting (keeping) the geometry as an attribute. [#586](https://github.com/rspatial/terra/issues/586) by Márcia Barbosa
- `saveRDS` and `serialize` methods for SpatRaster and SpatVector. [#549](https://github.com/rspatial/terra/issues/549) by Andrei Mîrț
- `xFromCol` and `yFromCol` now have a `<SpatRaster,missing>` method. [#583](https://github.com/rspatial/terra/issues/583) by Michael Sumner.
- `svc<sf>` method to deal with GeometryCollection types. [#585](https://github.com/rspatial/terra/issues/585) by Sarah Endicott
- `as.points<SpatRaster>` and `as.polygons<SpatRaster>` have a new argument `na.all=FALSE` that affects the interpretation of `na.rm`. [#548](https://github.com/rspatial/terra/issues/548) by Jean-Luc Dupouey.
- `setGDALconfig` and `getGDALconfig` to set GDAL configuration options. [#608](https://github.com/rspatial/terra/issues/608) by Erik Bolch.
- new argument `circular` to `rapp` to allow the start to be after the end (for if layers represent days of the year)
- new method `costDistance<SpatRaster>` 
- new methods `where.min` and `where.max` for `SpatRaster` to get the cell numbers for the extreme values in a SpatRaster. 
- new method `emptyGeoms<SpatVector>` to get the indices of empty (null) geometries
- new method `rasterizeGeom` to rasterize the geometry count or the area of (small) polygons or the length of lines.
- new method `not.na` for `SpatRaster` which is a shortcut for `!is.na(x)`.
- `as.list` implemented for `<SpatRasterDataset>`.
- `sources` implemented for `<SpatRasterDataset>`, `<SpatVector>` and `<SpatVectorProxy>` [#638](https://github.com/rspatial/terra/issues/638) by Andrew Gene Brown

## name changes

- delauny -> delaunay [#627](https://github.com/rspatial/terra/issues/627) by Derek Friend



# version 1.5-21

Released on 2022-02-17

- `writeVector` and `vect` now work with GPGK if the path has non-ascii characters [#518](https://github.com/rspatial/terra/issues/518)
- The results of `predict` with `cores > 1` and more than one output variable were garbled
- `zonal` dropped category names when using an external (R) function [#527](https://github.com/rspatial/terra/issues/527) by Jakub Nowosad
- focal/focalCpp showed strange patterns when the window size was larger than the block size [#519](https://github.com/rspatial/terra/issues/519) by Alex Ilich
- using `xy=TRUE` in `as.data.frame` normalized the names [#538](https://github.com/rspatial/terra/issues/538) by Kodi Arfer
- new argument `options` to `vrt` [#629](https://github.com/rspatial/terra/issues/629) by Monika Tomaszewska.


## enhancements 

- `makeTiles` has new arguments `extend` and `na.rm` [#520](https://github.com/rspatial/terra/issues/520) by by L. Dalby
- `project<SpatRaster>` now uses nearest neighbor as default method for RGB rasters
- new argument `na.rm=TRUE` to `unique`. [#561](https://github.com/rspatial/terra/issues/561) by Matthieu Stigler


# version 1.5-17

Released on 2022-01-30

## bug fixes

- `app<SpatRasterDataset>` ignored the filename. [#498](https://github.com/rspatial/terra/issues/498) by jszhao
- `vect<data.frame>` failed silently if xy coordinates were integers [#496](https://github.com/rspatial/terra/issues/496) by Márcia Barbosa
- The output of `aggregate<SpatRaster>` was malformed when `nrow(x) %% fact != 0`. [#492](https://github.com/rspatial/terra/issues/492) by Jean-François Bourdon
- Integer `NA`s in SpatVector attributes where only recognized on Windows [#491](https://github.com/rspatial/terra/issues/491) by Márcia Barbosa
- `plot<SpatVector>` failed when using a character variable with many unique values. [#489](https://github.com/rspatial/terra/issues/489) by Márcia Barbosa
- `rotate` failed on large files. Reported by Ujjawal Singh
- writing raster files with a color table could lead to a crash [#501](https://github.com/rspatial/terra/issues/501) by Kodi Arfer
- `crds` replicated the coordinates [#504](https://github.com/rspatial/terra/issues/504) by Murray Efford
- `as.data.frame<SpatRaster>` returned integers if the file stored values as integers, even if there was a scale/offset that creates decimal numbers [#509](https://github.com/rspatial/terra/issues/509) by Kodi Arfer
- `project` opened the input raster file in read/write mode intead of read mode. That did not work with files that cannot be updated.

## enhancements 

- `distance`, `gridDistance`, `direction` and `patches` now process all layers of the input SpatRaster. [#503](https://github.com/rspatial/terra/issues/503) by Chris Haak
- consistent copy-on-modify behavior in `()<-` methods. in-place updating available with `set.` methods such as `set.names` and `set.values`. [#493](https://github.com/rspatial/terra/issues/493) by Jean Romain and [511](https://github.com/rspatial/terra/issues/511) by Bryan Fuentes
- much faster writing of GPGK vector data by using a single transaction (following sf) [#460](https://github.com/rspatial/terra/issues/489) by Krzysztof Dyba
- `aggregate<SpatRaster>` now accepts functions that return more than one value per aggregated cell
- `writeVector` has new argument `insert` to add a layer to an existing file (e.g. GPKG).


## new

- new option `method="weights"` for `spatSample<SpatRaster>`
- new `mask<SpatVector,SpatVector>` method to select intersecting geometries
- new method `is.related`
- `values<SpatRaster>` has new option `na.rm=TRUE`. [#490](https://github.com/rspatial/terra/issues/490) by Henk Harmsen
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

Released on 2022-01-13

## bug fixes

- `setValues` and `init` failed (or even crashed R) when using a single value on a largish raster. [#414](https://github.com/rspatial/terra/issues/414)
- conversion from `sfc` to `SpatVector` lost the crs. [#415](https://github.com/rspatial/terra/issues/415) by Jean-Luc Dupouey
- `buffer` on a SpatRaster with no values caused a crash [#416](https://github.com/rspatial/terra/issues/416) by Sebastian Brinkmann
- `writeVector` now assumes "traditional GIS order" (long/lat) if the CRS specifies lat/long. [#333](
https://github.com/rspatial/terra/issues/333) by Agustin Lobo
- argument `main` was ignored in `density` when using a single layer SpatRaster [#424](https://github.com/rspatial/terra/issues/424) by dvictori
- Summary type math functions such as `min` and `mean`, when used with multiple SpatRasters and numbers, ignored additional SpatRasters [#426](https://github.com/rspatial/terra/issues/426) by Zhuonan Wang
- names are now conserved when creating a SpatRaster from a RasterStack that points to file(s) [#430](https://github.com/rspatial/terra/issues/430) by Dan Baston
- `classify` with `right=FALSE` ignored `include.lowest=TRUE` [#442](https://github.com/rspatial/terra/issues/442) by Alex Ilich
- `patches` now combines patches that connect across the data line [#366](https://github.com/rspatial/terra/issues/366) by Hirscht
- `patches(directions=8)` now connects in NE/SW direction [#451](https://github.com/rspatial/terra/issues/451) by Jean-François Bourdon.
- `centroids` now considers cases where SpatVector parts are nearest to each other when crossing the date line in stead of the zero-meridian [#366](https://github.com/rspatial/terra/issues/366) by Hirscht 
- `terrain` created empty (`NA`) rows between chunks used for processing large rasters. [#453](https://github.com/rspatial/terra/issues/452) by Robert Ritson.
- `inset` did not draw the "box" correctly. [#457](https://github.com/rspatial/terra/issues/457) by Márcia Barbosa
- `as.lines` now works with a points SpatVector [#464](https://github.com/rspatial/terra/issues/464) by Márcia Barbosa 


## enhancements 

- `values(x)<-` now accepts (hex coded) colors as values
- `focal` now wraps around the dateline like raster::focal [#242](https://github.com/rspatial/terra/issues/242) by Alexander Marbler
- `aggregate` now does not show a progress bar in all cases [#249](https://github.com/rspatial/terra/issues/249) by Lachlan
- `as.data.frame<SpatRaster> or <SpatVector>` are now also implemented as S3 methods to assure correct dispatch by other S3 methods such as `data.table::as.data.table`. See [#284](https://github.com/rspatial/terra/issues/284) by Patrick Schratz
- `crs` now shows the correct authority if it is not EPSG. [#419](https://github.com/rspatial/terra/issues/419) by Matthew Williamson
- It now possible to add a SpatRaster to an empty SpatRaster (with no values), even if it has a different geometry, ignoring the empty SpatRaster [#421](https://github.com/rspatial/terra/issues/421) by Alex Ilich.
- `rast<filename>` has a new argument `lyrs` to subset the layers and open the file in one step.
- `rast<array>` now has a crs and extent argument. [#439](https://github.com/rspatial/terra/issues/439) by RS-eco
- `type="xyz"` is now default in `rast<data.frame>`. [#438](https://github.com/rspatial/terra/issues/438) by RS-eco
- `classify` has a new argument `brackets` to show if a side of an interval is open or closed.
- further support for categorical data in `freq` and `as.data.frame`. [#441](https://github.com/rspatial/terra/issues/441) ngould7
- speed up in processing of multi-layer in memory data. [#437](https://github.com/rspatial/terra/issues/437) by Krzysztof Dyba
- `vect<matrix>` and `vect<data.frame>` are now much faster. [#413](https://github.com/rspatial/terra/issues/413) by BastienFR 	
- `extract` with points provided as a matrix or cell numbers is not much faster. [#341](https://github.com/rspatial/terra/issues/341)
- `focal` has a new argument `na.policy` that can be set to one of "all" (default), "only" or "omit". argument `na.only` has been removed, as you can now use `na.policy="only"`
- `inset` argument `border` changed to `perimeter` to allow passing `border` on to `plot<Spat*>`. [#456](https://github.com/rspatial/terra/issues/456) by Márcia Barbosa
- The compile-time and run-time versions of GEOS are now compared and a warning is given if they are not the same. [#459](https://github.com/rspatial/terra/issues/459) by Edzer Pebesma
- it is now possible to add sub-datasets to GPKG and GTiff files. [#300](https://github.com/rspatial/terra/issues/300) by gtitov
- general option `memfrac` can now be set to zero (in stead of not lower than 0.1). [#476](https://github.com/rspatial/terra/issues/476) by Matt Strimas-Mackey
- new argument `allowGaps` in `patches` to disallow gaps between patch IDs. See [#478](https://github.com/rspatial/terra/issues/478) by Dunbar Carpenter.


## new 

- timestamps and units are now saved to an auxiliary file (filename.aux.json) for all raster formats except NetCDF when using writeCDF (because in that case they are stored in the netcdf file)
- new method `mergeTime` to combine multiple rasters, perhaps partly overlapping in time, into a single time series
- new method `fillTime` that can add empty layers in between existing layers to assure that the time step between layers is constant 
- new method `approximate` to fill in missing values by cell across layers
- new methods `is.bool` and `as.bool` for SpatRaster and explicit recognition of Boolean raster data in various places (e.g., extract, plot)
- new methods `is.int` and `as.int` for SpatRaster. 
- when assigning integer values to a SpatRaster, or when reading an integer file, the corresponding layers are now classified as being of integer type [#446](https://github.com/rspatial/terra/issues/446) by L. Dalby
- new method `layerCor` (like `raster::layerStats`). [#420](https://github.com/rspatial/terra/issues/420) by Alex Ilich
- new method `focalCor` (like `raster::corLocal`). [#427](https://github.com/rspatial/terra/issues/427) by Zhuonan Wang
- new method `all.equal` for `SpatRaster`. See [#428](https://github.com/rspatial/terra/issues/428) by Dongdong Kong
- new method `math` for `SpatRaster` that implements the Math-generic methods *and* accepts a filename
- new method `sds<array>` 
- new method `rasterize<matrix>`, see [#413](https://github.com/rspatial/terra/issues/413) by BastienFR 	
- new method `colorize` to transform color representations 	
- new method `arrow` to draw a (North) arrow on a map. [#461](https://github.com/rspatial/terra/issues/461) by Márcia Barbosa
- new method `densify` to insert nodes between existing nodes of a line or polygon SpatVector
- new method `direction` for SpatRaster. [#462](https://github.com/rspatial/terra/issues/462) by Márcia Barbosa
- new method `focal3D` to compute focal values for a three-dimensional (row, column, layer) window
- new function `makeVRT` to create a vrt file for a file that needs a header to be read.
- new option `method="stratified"` for `spatSample<SpatRaster>`. [#470](https://github.com/rspatial/terra/issues/470) by Michael Mahoney
- new general option `memmax` to cap the amount of RAM that terra can be used in raster processing [#476](https://github.com/rspatial/terra/issues/476) by Matt Strimas-Mackey
- new method `gridDistance` to compute distances traversing a raster, perhaps with obstacles. [#477](https://github.com/rspatial/terra/issues/477) by Márcia Barbosa


# version 1.4-22

Released on 2021-11-24

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

Released on 2021-11-16

## bug fixes

- `terra` did not install with versions of GDAL below 3 [#402] by Alex Ilich.
- `distance` between two SpatVectors or matrices with `pairwise=FALSE` returned a matrix that was filled by column instead of by row [#403] by Paul Smith


# version 1.4-19

Released on 2021-11-15

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

Released on 2021-10-11

## enhancements

- the definition of `setValues` now has two arguments (`x` and `values`), just like `raster` had; to avoid reverse dependency problems with `raster`

# version 1.4-9

Released on 2021-10-07

## name changes

To avoid name conflicts with `sp` (via `raster`) `disaggregate` is now called `disagg` and `bbox,SpatRaster` and `bbox<SpatVector>` have been removed (but could be resurrected in `raster` or under another name).

## enhancements

- `project` and `resample` now choose the resampling method based on the first layer, using "near" for categorical data. Thanks to Matthew Lewis [#355]

## bug fixes

- `hist` failed with small samples. Issue [#356](https://github.com/rspatial/terra/issues/356) by Martin Queinnec


# version 1.4-7

Released on 2021-10-05

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

Released on 2021-08-20

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

Released on 2021-06-20

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

Released on 2021-05-13

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

Released on 2021-04-30

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

Released on 2021-04-14

## major changes 

- `c<SpatVector>` now returns a list. `rbind` is used to append SpatVector objects
- overhaul of handling of factors. `rats` has been removed, and `levels` and `cats` have changed


# version 1.1-4

- No news recorded for this version and earlier versions

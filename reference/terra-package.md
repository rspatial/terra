# Description of the methods in the terra package

`terra` provides methods to manipulate geographic (spatial) data in
"raster" and "vector" form. Raster data divide space into rectangular
grid cells and they are commonly used to represent spatially continuous
phenomena, such as elevation or the weather. Satellite images also have
this data structure, and in that context grid cells are often referred
to as pixels. In contrast, "vector" spatial data (points, lines,
polygons) are typically used to represent discrete spatial entities,
such as a road, country, or bus stop.

The package implements two main classes (data types): `SpatRaster` and
`SpatVector`. `SpatRaster` supports handling large raster files that
cannot be loaded into memory; local, focal, zonal, and global raster
operations; polygon, line and point to raster conversion; integration
with modeling methods to make spatial predictions; and more.
`SpatVector` supports all types of geometric operations such as
intersections.

Additional classes include `SpatExtent`, which is used to define a
spatial extent (bounding box); `SpatRasterDataset`, which represents a
collection of sub-datasets for the same area. Each sub-dataset is a
SpatRaster with possibly many layers, and may, for example, represent
different weather variables; and `SpatRasterCollection` and
`SpatVectorCollection` that are equivalent to lists of `SpatRaster` or
`SpatVector` objects. There is also a `SpatGraticule` class to assist in
adding a longitude/latitude lines and labels to a map with another
coordinate reference system.

These classes hold a C++ pointer to the data "reference class" and that
creates some limitations. They cannot be recovered from a saved R
session either or directly passed to nodes on a computer cluster.
Generally, you should use
[`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)
to save `SpatRaster` objects to disk (and pass a filename or cell values
to cluster nodes). Also see
[`wrap`](https://rspatial.github.io/terra/reference/wrap.md) and
[`saveRDS`](https://rspatial.github.io/terra/reference/serialize.md).
You should not write scripts that directly access this pointer, as its
user-interface is not stable.

The "terra" package is a replacement of the "raster" package. "terra"
has a very similar, but simpler, interface; it is faster, and it can do
much more. At the bottom of this page there is a table that shows
differences in the methods between the two packages.

Below is a list of some of the most important methods grouped by theme.

———————————————————————————————————————

## **SpatRaster**

———————————————————————————————————————

## I. Creating, combining and sub-setting

|                                                                                  |                                                               |
|----------------------------------------------------------------------------------|---------------------------------------------------------------|
| [`rast`](https://rspatial.github.io/terra/reference/rast.md)                     | Create a SpatRaster from scratch, file, or another object     |
| [`c`](https://rspatial.github.io/terra/reference/c.md)                           | Combine SpatRasters (multiple layers)                         |
| `add<-`                                                                          | Add a SpatRaster to another one                               |
| [`subset`](https://rspatial.github.io/terra/reference/subset.md) or `[[`, or `$` | Select layers of a SpatRaster                                 |
| [`selectRange`](https://rspatial.github.io/terra/reference/selectRange.md)       | Select cell values from different layers using an index layer |
| —————————                                                                        | ——————————————————————————————                                |

## II. Changing the spatial extent or resolution

Also see the methods in section VIII

|                                                                        |                                                                                    |
|------------------------------------------------------------------------|------------------------------------------------------------------------------------|
| [`merge`](https://rspatial.github.io/terra/reference/merge.md)         | Combine SpatRasters with different extents (but same origin and resolution)        |
| [`mosaic`](https://rspatial.github.io/terra/reference/mosaic.md)       | Combine SpatRasters with different extents using a function for overlapping cells  |
| [`crop`](https://rspatial.github.io/terra/reference/crop.md)           | Select a geographic subset of a SpatRaster                                         |
| [`extend`](https://rspatial.github.io/terra/reference/extend.md)       | Add rows and/or columns to a SpatRaster                                            |
| [`trim`](https://rspatial.github.io/terra/reference/trim.md)           | Trim a SpatRaster by removing exterior rows and/or columns that only have NAs      |
| [`aggregate`](https://rspatial.github.io/terra/reference/aggregate.md) | Combine cells of a SpatRaster to create larger cells                               |
| [`disagg`](https://rspatial.github.io/terra/reference/disaggregate.md) | Subdivide cells                                                                    |
| [`resample`](https://rspatial.github.io/terra/reference/resample.md)   | Resample (warp) values to a SpatRaster with a different origin and/or resolution   |
| [`project`](https://rspatial.github.io/terra/reference/project.md)     | Project (warp) values to a SpatRaster with a different coordinate reference system |
| [`shift`](https://rspatial.github.io/terra/reference/shift.md)         | Adjust the location of SpatRaster                                                  |
| [`flip`](https://rspatial.github.io/terra/reference/flip.md)           | Flip values horizontally or vertically                                             |
| [`rotate`](https://rspatial.github.io/terra/reference/rotate.md)       | Rotate values around the date-line (for lon/lat data)                              |
| [`t`](https://rspatial.github.io/terra/reference/transpose.md)         | Transpose a SpatRaster                                                             |
| —————————                                                              | ——————————————————————————————                                                     |

## III. Local (cell based) methods

### Apply-like methods

|                                                              |                                                                                                                                                                    |
|--------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`app`](https://rspatial.github.io/terra/reference/app.md)   | Apply a function to all cells, across layers, typically to summarize (as in [`base::apply`](https://rdrr.io/r/base/apply.html))                                    |
| [`tapp`](https://rspatial.github.io/terra/reference/tapp.md) | Apply a function to groups of layers (as in [`base::tapply`](https://rdrr.io/r/base/tapply.html) and [`stats::aggregate`](https://rdrr.io/r/stats/aggregate.html)) |
| [`lapp`](https://rspatial.github.io/terra/reference/lapp.md) | Apply a function to using the layers of a SpatRaster as variables                                                                                                  |
| [`sapp`](https://rspatial.github.io/terra/reference/sapp.md) | Apply a function to each layer                                                                                                                                     |
| [`rapp`](https://rspatial.github.io/terra/reference/rapp.md) | Apply a function to a spatially variable range of layers                                                                                                           |
| —————————                                                    | ——————————————————————————————                                                                                                                                     |

### Arithmetic, logical, and standard math methods

|                                                                                       |                                                                              |
|---------------------------------------------------------------------------------------|------------------------------------------------------------------------------|
| [`Arith-methods`](https://rspatial.github.io/terra/reference/arith-generic.md)        | Standard arithmetic methods (`+, -, *, ^, %%, %/%, /`)                       |
| [`Compare-methods`](https://rspatial.github.io/terra/reference/compare-generics.md)   | Comparison methods for SpatRaster (`==, !=, >, <, <=, >=m is.na, is.finite`) |
| [`not.na`](https://rspatial.github.io/terra/reference/not.na.md)                      | a one-step equivalent to `!is.na`                                            |
| [`Summary-methods`](https://rspatial.github.io/terra/reference/summarize-generics.md) | `mean, max, min, median, sum, range, prod,`                                  |
|                                                                                       | `any, all, stdev, which.min, which.max, anyNA, noNA, allNA`                  |
| [`Logic-methods`](https://rspatial.github.io/terra/reference/compare-generics.md)     | Boolean methods (`!, &, |`)                                                  |
| [`Math-methods`](https://rspatial.github.io/terra/reference/math-generics.md)         | `abs, sign, sqrt, ceiling, floor, trunc, cummax, cummin, cumprod,`           |
|                                                                                       | `cumsum, log, log10, log2, log1p, acos, acosh, asin, asinh, atan, atanh,`    |
|                                                                                       | `exp, expm1, cos, cosh, sin, sinh, tan, tanh, round, signif`                 |
| [`as.bool`](https://rspatial.github.io/terra/reference/is.bool.md)                    | create a Boolean (logical) SpatRaster                                        |
| [`as.int`](https://rspatial.github.io/terra/reference/is.bool.md)                     | create an integer (whole numbers) SpatRaster                                 |
| —————————                                                                             | ——————————————————————————————                                               |

### Other methods

|                                                                            |                                                                      |
|----------------------------------------------------------------------------|----------------------------------------------------------------------|
| [`approximate`](https://rspatial.github.io/terra/reference/approximate.md) | Compute missing values for cells by interpolation across layers      |
| [`roll`](https://rspatial.github.io/terra/reference/roll.md)               | Rolling functions such as the rolling mean                           |
| [`clamp`](https://rspatial.github.io/terra/reference/clamp.md)             | Restrict cell values to a minimum and/or maximum value               |
| [`cellSize`](https://rspatial.github.io/terra/reference/cellSize.md)       | Compute the area of cells                                            |
| [`classify`](https://rspatial.github.io/terra/reference/classify.md)       | (Re-)classify values                                                 |
| [`subst`](https://rspatial.github.io/terra/reference/subst.md)             | Substitute (replace) cell values                                     |
| [`cover`](https://rspatial.github.io/terra/reference/cover.md)             | First layer covers second layer except where the first layer is `NA` |
| [`init`](https://rspatial.github.io/terra/reference/init.md)               | Initialize cells with new values                                     |
| [`mask`](https://rspatial.github.io/terra/reference/mask.md)               | Replace values in a SpatRaster based on values in another SpatRaster |
| [`which.lyr`](https://rspatial.github.io/terra/reference/which.md)         | which is the first layer that is `TRUE`?                             |
| [`segregate`](https://rspatial.github.io/terra/reference/segregate.md)     | Make a 0/1 layer for each unique value                               |
| [`rangeFill`](https://rspatial.github.io/terra/reference/rangeFill.md)     | Make a 0/1 SpatRaster for a time series                              |
| [`regress`](https://rspatial.github.io/terra/reference/regress.md)         | Cell-based regression models                                         |
| —————————                                                                  | ——————————————————————————————                                       |

## IV. Zonal and global methods

|                                                                      |                                                            |
|----------------------------------------------------------------------|------------------------------------------------------------|
| [`expanse`](https://rspatial.github.io/terra/reference/expanse.md)   | Compute the summed area of cells                           |
| [`crosstab`](https://rspatial.github.io/terra/reference/crosstab.md) | Cross-tabulate two SpatRasters                             |
| [`freq`](https://rspatial.github.io/terra/reference/freq.md)         | Frequency table of SpatRaster cell values                  |
| [`global`](https://rspatial.github.io/terra/reference/global.md)     | Summarize SpatRaster cell values with a function           |
| [`quantile`](https://rspatial.github.io/terra/reference/quantile.md) | Quantiles                                                  |
| [`layerCor`](https://rspatial.github.io/terra/reference/layerCor.md) | Correlation between layers                                 |
| [`stretch`](https://rspatial.github.io/terra/reference/stretch.md)   | Stretch values                                             |
| [`scale`](https://rspatial.github.io/terra/reference/scale.md)       | Scale values                                               |
| [`summary`](https://rspatial.github.io/terra/reference/summary.md)   | Summary of the values of a SpatRaster (quartiles and mean) |
| [`unique`](https://rspatial.github.io/terra/reference/unique.md)     | Get the unique values in a SpatRaster                      |
| [`zonal`](https://rspatial.github.io/terra/reference/zonal.md)       | Summarize a SpatRaster by zones in another SpatRaster      |
| —————————                                                            | ——————————————————————————————                             |

## V. Situation (spatial context) based methods

|                                                                          |                                                                                       |
|--------------------------------------------------------------------------|---------------------------------------------------------------------------------------|
| [`adjacent`](https://rspatial.github.io/terra/reference/adjacent.md)     | Identify cells that are adjacent to a set of cells of a SpatRaster                    |
| [`boundaries`](https://rspatial.github.io/terra/reference/boundaries.md) | Detection of boundaries (edges)                                                       |
| [`distance`](https://rspatial.github.io/terra/reference/distance.md)     | Shortest distance to a cell that is not `NA` or to or from a vector object            |
| [`gridDist`](https://rspatial.github.io/terra/reference/gridDist.md)     | Shortest distance through adjacent grid cells                                         |
| [`costDist`](https://rspatial.github.io/terra/reference/costDist.md)     | Shortest distance considering cell-varying friction                                   |
| [`direction`](https://rspatial.github.io/terra/reference/direction.md)   | Direction (azimuth) to or from cells that are not `NA`                                |
| [`focal`](https://rspatial.github.io/terra/reference/focal.md)           | Focal (neighborhood; moving window) functions                                         |
| [`focal3D`](https://rspatial.github.io/terra/reference/focal3D.md)       | Three dimensional (row, col, lyr) focal functions                                     |
| [`focalCpp`](https://rspatial.github.io/terra/reference/focalCpp.md)     | Faster focal by using custom C++ functions                                            |
| [`focalReg`](https://rspatial.github.io/terra/reference/focalReg.md)     | Regression between layers for focal areas                                             |
| [`focalPairs`](https://rspatial.github.io/terra/reference/focalPairs.md) | Apply a function (e.g. a correlation coefficient) to focal values for pairs of layers |
| [`patches`](https://rspatial.github.io/terra/reference/patches.md)       | Find patches (clumps)                                                                 |
| [`sieve`](https://rspatial.github.io/terra/reference/sieve.md)           | Sieve filter to remove small patches                                                  |
| [`terrain`](https://rspatial.github.io/terra/reference/terrain.md)       | Compute slope, aspect and other terrain characteristics from elevation data           |
| [`viewshed`](https://rspatial.github.io/terra/reference/viewshed.md)     | Compute viewshed (showing areas that are visible from a particular location           |
| [`shade`](https://rspatial.github.io/terra/reference/shade.md)           | Compute hill shade from slope and aspect layers                                       |
| [`autocor`](https://rspatial.github.io/terra/reference/autocor.md)       | Compute global or local spatial autocorrelation                                       |
| —————————                                                                | ——————————————————————————————                                                        |

## VI. Model predictions

|                                                                                                                                             |                                                                            |
|---------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------|
| [`predict`](https://rspatial.github.io/terra/reference/predict.md)                                                                          | Predict a non-spatial (regression or classification) model to a SpatRaster |
| [`interpolate`](https://rspatial.github.io/terra/reference/interpolate.md)                                                                  | Predict a spatial model to a SpatRaster                                    |
| [`interpIDW`](https://rspatial.github.io/terra/reference/interpIDW.md)                                                                      | Inverse-distance-weighted interpolation                                    |
| [`interpNear`](https://rspatial.github.io/terra/reference/interpNear.md)                                                                    | Nearest neighbor interpolation                                             |
| [`k_means`](https://rspatial.github.io/terra/reference/k_means.md)                                                                          | k-means clustering of SpatRaster data                                      |
| [`princomp`](https://rspatial.github.io/terra/reference/princomp.md)` and `[`prcomp`](https://rspatial.github.io/terra/reference/prcomp.md) | Principal Component Analysis (PCA) with raster data                        |
| —————————                                                                                                                                   | ——————————————————————————————                                             |

## VII. Accessing cell values

Apart from the function listed below, you can also use indexing with `[`
with cell numbers, and row and/or column numbers  

|                                                                                |                                                                                                    |
|--------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------|
| [`values`](https://rspatial.github.io/terra/reference/values.md)               | cell values (fails with very large rasters)                                                        |
| `values<-`                                                                     | Set new values to the cells of a SpatRaster                                                        |
| [`setValues`](https://rspatial.github.io/terra/reference/setValues.md)         | Set new values to the cells of a SpatRaster                                                        |
| [`as.matrix`](https://rspatial.github.io/terra/reference/coerce.md)            | Get cell values as a matrix                                                                        |
| [`as.array`](https://rspatial.github.io/terra/reference/coerce.md)             | Get cell values as an array                                                                        |
| [`as.data.frame`](https://rspatial.github.io/terra/reference/as.data.frame.md) | get cell values as a data.frame (including class lables)                                           |
| [`extract`](https://rspatial.github.io/terra/reference/extract.md)             | Extract cell values from a SpatRaster (with cell numbers, coordinates, points, lines, or polygons) |
| [`extractAlong`](https://rspatial.github.io/terra/reference/extractAlong.md)   | Extract cell values along a line such that the values are in the right order                       |
| [`spatSample`](https://rspatial.github.io/terra/reference/sample.md)           | Take a sample (regular, random, stratified, weighted) sample from a SpatRaster                     |
| [`minmax`](https://rspatial.github.io/terra/reference/minmax.md)               | Get the minimum and maximum value of the cells of a SpatRaster (if known)                          |
| [`setMinMax`](https://rspatial.github.io/terra/reference/minmax.md)            | Compute the minimum and maximum value of a SpatRaster if these are not known                       |
| —————————                                                                      | ——————————————————————————————                                                                     |

## VIII. Getting and setting dimensions

Get or set basic parameters of SpatRasters. If there are values
associated with a SpatRaster (either in memory or via a link to a file)
these are lost when you change the number of columns or rows or the
resolution. This is not the case when the extent is changed (as the
number of columns and rows will not be affected). Similarly, with
**crs** you can set the coordinate reference system, but this does not
transform the data (see
[project](https://rspatial.github.io/terra/reference/project.md) for
that).

|                                                                            |                                                                                 |
|----------------------------------------------------------------------------|---------------------------------------------------------------------------------|
| [`ncol`](https://rspatial.github.io/terra/reference/dimensions.md)         | The number of columns                                                           |
| [`nrow`](https://rspatial.github.io/terra/reference/dimensions.md)         | The number of rows                                                              |
| [`ncell`](https://rspatial.github.io/terra/reference/dimensions.md)        | The number of cells (can not be set directly, only via ncol or nrow)            |
| [`res`](https://rspatial.github.io/terra/reference/dimensions.md)          | The resolution (x and y)                                                        |
| [`nlyr`](https://rspatial.github.io/terra/reference/dimensions.md)         | Get or set the number of layers                                                 |
| [`names`](https://rspatial.github.io/terra/reference/names.md)             | Get or set the layer names                                                      |
| [`xres`](https://rspatial.github.io/terra/reference/dimensions.md)         | The x resolution (can be set with res)                                          |
| [`yres`](https://rspatial.github.io/terra/reference/dimensions.md)         | The y resolution (can be set with res)                                          |
| [`xmin`](https://rspatial.github.io/terra/reference/xmin.md)               | The minimum x coordinate (or longitude)                                         |
| [`xmax`](https://rspatial.github.io/terra/reference/xmin.md)               | The maximum x coordinate (or longitude)                                         |
| [`ymin`](https://rspatial.github.io/terra/reference/xmin.md)               | The minimum y coordinate (or latitude)                                          |
| [`ymax`](https://rspatial.github.io/terra/reference/xmin.md)               | The maximum y coordinate (or latitude)                                          |
| [`ext`](https://rspatial.github.io/terra/reference/ext.md)                 | Get or set the extent (minimum and maximum x and y coordinates ("bounding box") |
| [`origin`](https://rspatial.github.io/terra/reference/origin.md)           | The origin of a SpatRaster                                                      |
| [`sources`](https://rspatial.github.io/terra/reference/sources.md)         | Get the filename(s) to which a SpatRaster is linked                             |
| [`inMemory`](https://rspatial.github.io/terra/reference/sources.md)        | Are the data sources in memory (or on disk)?                                    |
| [`toMemory`](https://rspatial.github.io/terra/reference/toMemory.md)       | Force data sources to memory (not recommended)?                                 |
| [`compareGeom`](https://rspatial.github.io/terra/reference/compareGeom.md) | Compare the geometry of SpatRasters                                             |
| [`NAflag`](https://rspatial.github.io/terra/reference/NAflag.md)           | Set the `NA` value (for reading from a file with insufficient metadata)         |
| —————————                                                                  | ——————————————————————————————                                                  |

## IX. Computing row, column, cell numbers and coordinates

Cell numbers start at 1 in the upper-left corner. They increase within
rows, from left to right, and then row by row from top to bottom.
Likewise, row numbers start at 1 at the top of the raster, and column
numbers start at 1 at the left side of the raster.

|                                                                                     |                                                              |
|-------------------------------------------------------------------------------------|--------------------------------------------------------------|
| [`xFromCol`](https://rspatial.github.io/terra/reference/xyCellFrom.md)              | x-coordinates from column numbers                            |
| [`yFromRow`](https://rspatial.github.io/terra/reference/xyCellFrom.md)              | y-coordinates from row numbers                               |
| [`xFromCell`](https://rspatial.github.io/terra/reference/xyCellFrom.md)             | x-coordinates from row numbers                               |
| [`yFromCell`](https://rspatial.github.io/terra/reference/xyCellFrom.md)             | y-coordinates from cell numbers                              |
| [`xyFromCell`](https://rspatial.github.io/terra/reference/xyCellFrom.md)            | x and y coordinates from cell numbers                        |
| [`colFromX`](https://rspatial.github.io/terra/reference/xyCellFrom.md)              | Column numbers from x-coordinates (or longitude)             |
| [`rowFromY`](https://rspatial.github.io/terra/reference/xyCellFrom.md)              | Row numbers from y-coordinates (or latitude)                 |
| [`rowColFromCell`](https://rspatial.github.io/terra/reference/xyCellFrom.md)        | Row and column numbers from cell numbers                     |
| [`cellFromXY`](https://rspatial.github.io/terra/reference/xyCellFrom.md)            | Cell numbers from x and y coordinates                        |
| [`cellFromRowCol`](https://rspatial.github.io/terra/reference/xyCellFrom.md)        | Cell numbers from row and column numbers                     |
| [`cellFromRowColCombine`](https://rspatial.github.io/terra/reference/xyCellFrom.md) | Cell numbers from all combinations of row and column numbers |
| [`cells`](https://rspatial.github.io/terra/reference/cells.md)                      | Cell numbers for a SpatVector or SpatExtent                  |
| —————————                                                                           | ——————————————————————————————                               |

## X. Depth related methods

`depth` can be used to explicitly a third or fourth dimension of a
SpatRaster.

|                                                                    |                                      |
|--------------------------------------------------------------------|--------------------------------------|
| [`depth`](https://rspatial.github.io/terra/reference/depth.md)     | Get or set depth dimension values () |
| [`depthName`](https://rspatial.github.io/terra/reference/depth.md) | Set or get the depth name            |
| [`depthUnit`](https://rspatial.github.io/terra/reference/depth.md) | Set or get the depth unit            |
| —————————                                                          | ——————————————————————————————       |

## XI. Time related methods

`time` can be used to explicitly a third or fourth dimension of a
SpatRaster.

|                                                                        |                                                                                                         |
|------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
| [`time`](https://rspatial.github.io/terra/reference/time.md)           | Get or set time                                                                                         |
| [`fillTime`](https://rspatial.github.io/terra/reference/fillTime.md)   | can add empty layers in between existing layers to assure that the time step between layers is constant |
| [`mergeTime`](https://rspatial.github.io/terra/reference/mergeTime.md) | combine multiple rasters, perhaps partly overlapping in time, into a single time series                 |
| —————————                                                              | ——————————————————————————————                                                                          |

## XII. Methods for categorical rasters

|                                                                        |                                                                |
|------------------------------------------------------------------------|----------------------------------------------------------------|
| [`is.factor`](https://rspatial.github.io/terra/reference/is.bool.md)   | Are there categorical layers?                                  |
| [`levels`](https://rspatial.github.io/terra/reference/factors.md)      | Get active categories, or set categories                       |
| [`activeCat`](https://rspatial.github.io/terra/reference/activeCat.md) | Get or set the active category                                 |
| [`cats`](https://rspatial.github.io/terra/reference/factors.md)        | Get categories (active and inactive)                           |
| [`set.cats`](https://rspatial.github.io/terra/reference/inplace.md)    | Set categories in place                                        |
| [`concats`](https://rspatial.github.io/terra/reference/concats.md)     | Combine SpatRasters with different categories                  |
| [`catalyze`](https://rspatial.github.io/terra/reference/catalyze.md)   | Create a layer for each category                               |
| [`as.numeric`](https://rspatial.github.io/terra/reference/catalyze.md) | use the active category to create a non-categorical SpatRaster |
| [`as.factor`](https://rspatial.github.io/terra/reference/is.bool.md)   | Make the layers of a SpatRaster categorical                    |
| —————————                                                              | ——————————————————————————————                                 |

## XIII. Writing SpatRaster files

### Basic

|                                                                            |                                                                                          |
|----------------------------------------------------------------------------|------------------------------------------------------------------------------------------|
| [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md) | Write all values of SpatRaster to disk. You can set the filetype, datatype, compression. |
| [`writeCDF`](https://rspatial.github.io/terra/reference/writeCDF.md)       | Write SpatRaster data to a netCDF file                                                   |
| —————————                                                                  | ——————————————————————————————                                                           |

### Advanced

|                                                                          |                                                         |
|--------------------------------------------------------------------------|---------------------------------------------------------|
| [`readStart`](https://rspatial.github.io/terra/reference/readwrite.md)   | Open file connections for efficient multi-chunk reading |
| [`readValues`](https://rspatial.github.io/terra/reference/readwrite.md)  | Read some values from an opened file                    |
| [`readStop`](https://rspatial.github.io/terra/reference/readwrite.md)    | Close file connections                                  |
| [`writeStart`](https://rspatial.github.io/terra/reference/readwrite.md)  | Open a file for writing                                 |
| [`writeValues`](https://rspatial.github.io/terra/reference/readwrite.md) | Write some values to an opened file                     |
| [`writeStop`](https://rspatial.github.io/terra/reference/readwrite.md)   | Close the file after writing                            |
| [`blocks`](https://rspatial.github.io/terra/reference/readwrite.md)      | Get blocksize for reading files (when not writing)      |
| —————————                                                                | ——————————————————————————————                          |

## XIV. Miscellaneous SpatRaster methods

|                                                                              |                                                                                          |
|------------------------------------------------------------------------------|------------------------------------------------------------------------------------------|
| [`terraOptions`](https://rspatial.github.io/terra/reference/terraOptions.md) | Show, set, or get session options, mostly to control memory use and to set write options |
| [`sources`](https://rspatial.github.io/terra/reference/sources.md)           | Show the data sources of a SpatRaster                                                    |
| [`tmpFiles`](https://rspatial.github.io/terra/reference/tmpFile.md)          | Show or remove temporary files                                                           |
| [`mem_info`](https://rspatial.github.io/terra/reference/mem.md)              | memory needs and availability                                                            |
| [`inMemory`](https://rspatial.github.io/terra/reference/sources.md)          | Are the cell values in memory?                                                           |
| —————————                                                                    | ——————————————————————————————                                                           |

## XV. SpatRasterDataset

A SpatRasterDataset contains SpatRasters that represent sub-datasets for
the same area. They all have the same extent and resolution.

|                                                                |                                                                                           |
|----------------------------------------------------------------|-------------------------------------------------------------------------------------------|
| [`sds`](https://rspatial.github.io/terra/reference/sds.md)     | Create a SpatRasterDataset from a file with subdatasets (ncdf or hdf) or from SpatRasters |
| `[` or `$`                                                     | Extract a SpatRaster                                                                      |
| [`names`](https://rspatial.github.io/terra/reference/names.md) | Get the names of the sub-datasets                                                         |
| —————————                                                      | ——————————————————————————————                                                            |

## XVI. SpatRasterCollections

A SpatRasterCollection is a vector of SpatRaster objects. Unlike for a
SpatRasterDataset, there the extent and resolution of the SpatRasters do
not need to match each other.

|                                                                      |                                                                                            |
|----------------------------------------------------------------------|--------------------------------------------------------------------------------------------|
| [`sprc`](https://rspatial.github.io/terra/reference/sprc.md)         | create a SpatRasterCollection from (a list of) SpatRasters                                 |
| [`length`](https://rspatial.github.io/terra/reference/dimensions.md) | how many SpatRasters does the SpatRasterCollection have?                                   |
| [`crop`](https://rspatial.github.io/terra/reference/crop.md)         | crop a SpatRasterCollection                                                                |
| [`impose`](https://rspatial.github.io/terra/reference/impose.md)     | force the members of SpatRasterCollection to the same geometry                             |
| [`merge`](https://rspatial.github.io/terra/reference/merge.md)       | merge the members of a SpatRasterCollection                                                |
| [`mosaic`](https://rspatial.github.io/terra/reference/mosaic.md)     | mosaic (merge with a function for overlapping areas) the members of a SpatRasterCollection |
| [`[`](https://rspatial.github.io/terra/reference/subset_single.md)   | extract a SpatRaster                                                                       |
| —————————                                                            | ——————————————————————————————                                                             |

## **SpatVector**

———————————————————————————————————————

## XVII. Create SpatVector objects

|                                                                                |                                                                                         |
|--------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------|
| [`vect`](https://rspatial.github.io/terra/reference/vect.md)                   | Create a SpatVector from a file (for example a "shapefile") or from another object      |
| [`vector_layers`](https://rspatial.github.io/terra/reference/vector_layers.md) | list or delete layers in a vector database such as GPGK                                 |
| `rbind`                                                                        | append SpatVectors of the same geometry type                                            |
| [`unique`](https://rspatial.github.io/terra/reference/unique.md)               | remove duplicates                                                                       |
| [`na.omit`](https://rspatial.github.io/terra/reference/na.omit.md)             | remove empty geometries and/or fields that are `NA`                                     |
| [`project`](https://rspatial.github.io/terra/reference/project.md)             | Project a SpatVector to a different coordinate reference system                         |
| [`writeVector`](https://rspatial.github.io/terra/reference/writeVector.md)     | Write SpatVector data to disk                                                           |
| [`centroids`](https://rspatial.github.io/terra/reference/centroids.md)         | Get the centroids of a SpatVector                                                       |
| [`voronoi`](https://rspatial.github.io/terra/reference/voronoi.md)             | Voronoi diagram                                                                         |
| [`delaunay`](https://rspatial.github.io/terra/reference/voronoi.md)            | Delaunay triangles                                                                      |
| [`hull`](https://rspatial.github.io/terra/reference/convhull.md)               | Compute a convex, circular, or rectangular hull around the (geometries of) a SpatVector |
| [`fillHoles`](https://rspatial.github.io/terra/reference/fill.md)              | Remove or extract holes from polygons                                                   |
| —————————                                                                      | ——————————————————————————————                                                          |

## XVIII. Properties of SpatVector objects

|                                                                            |                                                                                      |
|----------------------------------------------------------------------------|--------------------------------------------------------------------------------------|
| [`geom`](https://rspatial.github.io/terra/reference/geometry.md)           | returns the geometries as matrix or WKT                                              |
| [`crds`](https://rspatial.github.io/terra/reference/crds.md)               | returns the coordinates as a matrix                                                  |
| [`ncol`](https://rspatial.github.io/terra/reference/dimensions.md)         | The number of columns (of the attributes)                                            |
| [`nrow`](https://rspatial.github.io/terra/reference/dimensions.md)         | The number of rows (of the geometries and attributes)                                |
| [`names`](https://rspatial.github.io/terra/reference/names.md)             | Get or set the layer names                                                           |
| [`ext`](https://rspatial.github.io/terra/reference/ext.md)                 | Get the extent (minimum and maximum x and y coordinates ("bounding box")             |
| [`crs`](https://rspatial.github.io/terra/reference/crs.md)                 | The coordinate reference system (map projection)                                     |
| [`linearUnits`](https://rspatial.github.io/terra/reference/linearUnits.md) | returns the linear units of the crs (in meter)                                       |
| [`is.lonlat`](https://rspatial.github.io/terra/reference/is.lonlat.md)     | Test if an object has (or may have) a longitude/latitude coordinate reference system |
| —————————                                                                  | ——————————————————————————————                                                       |

## XIX. Geometric queries

|                                                                      |                                                                           |
|----------------------------------------------------------------------|---------------------------------------------------------------------------|
| [`adjacent`](https://rspatial.github.io/terra/reference/adjacent.md) | find adjacent polygons                                                    |
| [`expanse`](https://rspatial.github.io/terra/reference/expanse.md)   | computes the area covered by polygons                                     |
| [`nearby`](https://rspatial.github.io/terra/reference/nearby.md)     | find nearby geometries                                                    |
| [`nearest`](https://rspatial.github.io/terra/reference/nearby.md)    | find the nearest geometries                                               |
| [`relate`](https://rspatial.github.io/terra/reference/relate.md)     | geometric relationships such as "intersects", "overlaps", and "touches"   |
| [`perim`](https://rspatial.github.io/terra/reference/perim.md)       | computes the length of the perimeter of polygons, and the length of lines |
| —————————                                                            | ——————————————————————————————                                            |

## XX. Geometric operations

|                                                                                |                                                              |
|--------------------------------------------------------------------------------|--------------------------------------------------------------|
| [`erase`](https://rspatial.github.io/terra/reference/erase.md) or "-"          | erase (parts of) geometries                                  |
| [`intersect`](https://rspatial.github.io/terra/reference/intersect.md) or "\*" | intersect geometries                                         |
| [`union`](https://rspatial.github.io/terra/reference/union.md) or "+"          | Merge geometries                                             |
| [`cover`](https://rspatial.github.io/terra/reference/cover.md)                 | update polygons                                              |
| [`symdif`](https://rspatial.github.io/terra/reference/symdif.md)               | symmetrical difference of two polygons                       |
| [`aggregate`](https://rspatial.github.io/terra/reference/aggregate.md)         | dissolve smaller polygons into larger ones                   |
| [`buffer`](https://rspatial.github.io/terra/reference/buffer.md)               | buffer geometries                                            |
| [`disagg`](https://rspatial.github.io/terra/reference/disaggregate.md)         | split multi-geometries into separate geometries              |
| [`crop`](https://rspatial.github.io/terra/reference/crop.md)                   | clip geometries using a rectangle (SpatExtent) or SpatVector |
| —————————                                                                      | ——————————————————————————————                               |

## XXI. SpatVector attributes

We use the term "attributes" for the tabular data (data.frame)
associated with vector geometries.

|                                                                                |                                                                            |
|--------------------------------------------------------------------------------|----------------------------------------------------------------------------|
| [`extract`](https://rspatial.github.io/terra/reference/extract.md)             | spatial queries between SpatVector and SpatVector (e.g. point in polygons) |
| [`spatSample`](https://rspatial.github.io/terra/reference/sample.md)           | Take a regular or random point sample from polygons or lines               |
| [`sel`](https://rspatial.github.io/terra/reference/select.md)                  | select - interactively select geometries                                   |
| [`click`](https://rspatial.github.io/terra/reference/click.md)                 | identify attributes by clicking on a map                                   |
| [`merge`](https://rspatial.github.io/terra/reference/merge.md)                 | Join a table with a SpatVector                                             |
| [`as.data.frame`](https://rspatial.github.io/terra/reference/as.data.frame.md) | get attributes as a data.frame                                             |
| [`as.list`](https://rspatial.github.io/terra/reference/as.list.md)             | get attributes as a list                                                   |
| [`values`](https://rspatial.github.io/terra/reference/values.md)               | Get the attributes of a SpatVector                                         |
| `values<-`                                                                     | Set new attributes to the geometries of a SpatRaster                       |
| [`sort`](https://rspatial.github.io/terra/reference/sort.md)                   | sort SpatVector by the values in a field                                   |
| —————————                                                                      | ——————————————————————————————                                             |

## XXII. Change geometries (for display, experimentation)

|                                                                    |                                                                                                      |
|--------------------------------------------------------------------|------------------------------------------------------------------------------------------------------|
| [`shift`](https://rspatial.github.io/terra/reference/shift.md)     | change the position geometries by shifting their coordinates in horizontal and/or vertical direction |
| [`spin`](https://rspatial.github.io/terra/reference/spin.md)       | rotate geometries around an origin                                                                   |
| [`rescale`](https://rspatial.github.io/terra/reference/rescale.md) | shrink (or expand) geometries, for example to make an inset map                                      |
| [`flip`](https://rspatial.github.io/terra/reference/flip.md)       | flip geometries vertically or horizontally                                                           |
| [`t`](https://rspatial.github.io/terra/reference/transpose.md)     | transpose geometries (switch x and y)                                                                |
| —————————                                                          | ——————————————————————————————                                                                       |

## XXIII. Geometry properties and topology

|                                                                                            |                                                                                                 |
|--------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| [`width`](https://rspatial.github.io/terra/reference/width.md)                             | the minimum diameter of the geometries                                                          |
| [`clearance`](https://rspatial.github.io/terra/reference/width.md)                         | the minimum clearance of the geometries                                                         |
| [`sharedPaths`](https://rspatial.github.io/terra/reference/sharedPaths.md)                 | shared paths (arcs) between line or polygon geometries                                          |
| [`simplifyGeom`](https://rspatial.github.io/terra/reference/simplify.md)                   | simplify geometries                                                                             |
| [`gaps`](https://rspatial.github.io/terra/reference/gaps.md)                               | find gaps between polygon geometries                                                            |
| [`fillHoles`](https://rspatial.github.io/terra/reference/fill.md)                          | get or remove the polygon holes                                                                 |
| [`makeNodes`](https://rspatial.github.io/terra/reference/topology.md)                      | create nodes on lines                                                                           |
| [`mergeLines`](https://rspatial.github.io/terra/reference/topology.md)                     | connect lines to form polygons                                                                  |
| [`removeDupNodes`](https://rspatial.github.io/terra/reference/topology.md)                 | remove duplicate nodes in geometries and optionally rounds the coordinates                      |
| [`is.valid`](https://rspatial.github.io/terra/reference/is.valid.md)                       | check if geometries are valid                                                                   |
| [`makeValid`](https://rspatial.github.io/terra/reference/is.valid.md)                      | attempt to repair invalid geometries                                                            |
| [`snap`](https://rspatial.github.io/terra/reference/topology.md)                           | make boundaries of geometries identical if they are very close to each other                    |
| [`erase`](https://rspatial.github.io/terra/reference/erase.md)` (single argument)`         | remove parts of geometries that overlap                                                         |
| [`union`](https://rspatial.github.io/terra/reference/union.md)` (single argument)`         | create new polygons such that there are no overlapping polygons                                 |
| [`rotate`](https://rspatial.github.io/terra/reference/rotate.md)                           | rotate to (dis-) connect them across the date-line                                              |
| [`normalize.longitude`](https://rspatial.github.io/terra/reference/normalize.longitude.md) | move geometries that are outside of the -180 to 180 degrees range.                              |
| [`elongate`](https://rspatial.github.io/terra/reference/elongate.md)                       | make lines longer by extending both sides                                                       |
| [`combineGeoms`](https://rspatial.github.io/terra/reference/combineGeoms.md)               | combine geometries that overlap, share a border, or are within a minimum distance of each other |
| [`forceCCW`](https://rspatial.github.io/terra/reference/forceCCW.md)                       | force counter-clockwise polygon winding                                                         |
| —————————                                                                                  | ——————————————————————————————                                                                  |

## XXIV. SpatVectorCollections

A SpatVectorCollection is a vector of SpatVector objects.

|                                                                      |                                                                   |
|----------------------------------------------------------------------|-------------------------------------------------------------------|
| [`svc`](https://rspatial.github.io/terra/reference/svc.md)           | create a SpatVectorCollection from (a list of) SpatVector objects |
| [`length`](https://rspatial.github.io/terra/reference/dimensions.md) | how many SpatRasters does the SpatRasterCollection have?          |
| [`[`](https://rspatial.github.io/terra/reference/subset_single.md)   | extract a SpatVector                                              |
| —————————                                                            | ——————————————————————————————                                    |

## XXV. Coordinate reference system method

|                                                                            |                                                                                      |
|----------------------------------------------------------------------------|--------------------------------------------------------------------------------------|
| [`crs`](https://rspatial.github.io/terra/reference/crs.md)                 | Get or set the coordinate reference system (map projection) of a Spat\* object       |
| [`is.lonlat`](https://rspatial.github.io/terra/reference/is.lonlat.md)     | Test if an object has (or may have) a longitude/latitude coordinate reference system |
| [`linearUnits`](https://rspatial.github.io/terra/reference/linearUnits.md) | returns the linear units of the crs (in meter)                                       |
| —————————                                                                  | ——————————————————————————————                                                       |

## **Other classes**

———————————————————————————————————————

## XXVI. SpatExtent

|                                                                               |                                                                                                                           |
|-------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------|
| [`ext`](https://rspatial.github.io/terra/reference/ext.md)                    | Create a SpatExtent object. For example to [`crop`](https://rspatial.github.io/terra/reference/crop.md) a Spatial dataset |
| [`intersect`](https://rspatial.github.io/terra/reference/intersect.md)        | Intersect two SpatExtent objects, same as `-`                                                                             |
| [`union`](https://rspatial.github.io/terra/reference/union.md)                | Combine two SpatExtent objects, same as `+`                                                                               |
| [`Math-methods`](https://rspatial.github.io/terra/reference/math-generics.md) | round/floor/ceiling of a SpatExtent                                                                                       |
| [`align`](https://rspatial.github.io/terra/reference/align.md)                | Align a SpatExtent with a SpatRaster                                                                                      |
| [`draw`](https://rspatial.github.io/terra/reference/draw.md)                  | Create a SpatExtent by drawing it on top of a map (plot)                                                                  |
| —————————                                                                     | ——————————————————————————————                                                                                            |

## XXVII. SpatGraticule

|                                                                                       |                                |
|---------------------------------------------------------------------------------------|--------------------------------|
| [`graticule`](https://rspatial.github.io/terra/reference/graticule.md)                | Create a graticule             |
| [`crop`](https://rspatial.github.io/terra/reference/crop.md)                          | crop a graticule               |
| [`plot<SpatGraticule>`](https://rspatial.github.io/terra/reference/plot_graticule.md) | plot a graticule               |
| —————————                                                                             | —————————————————————————————— |

## **General methods**

———————————————————————————————————————

## XXVIII. Conversion between spatial data objects from different packages

You can coerce SpatRasters to Raster\* objects, after loading the
`raster` package, with `as(object, "Raster")`, or `raster(object)` or
`brick(object)` or `stack(object)`

|                                                                                |                                                                         |
|--------------------------------------------------------------------------------|-------------------------------------------------------------------------|
| [`rast`](https://rspatial.github.io/terra/reference/rast.md)                   | SpatRaster from matrix and other objects                                |
| [`vect`](https://rspatial.github.io/terra/reference/vect.md)                   | SpatVector from `sf` or `Spatial*` vector data                          |
| [`sf::st_as_sf`](https://r-spatial.github.io/sf/reference/st_as_sf.html)       | sf object from SpatVector                                               |
| [`rasterize`](https://rspatial.github.io/terra/reference/rasterize.md)         | Rasterizing points, lines or polygons                                   |
| [`rasterizeWin`](https://rspatial.github.io/terra/reference/rasterizeWin.md)   | Rasterize points with a moving window                                   |
| [`rasterizeGeom`](https://rspatial.github.io/terra/reference/rasterizeGeom.md) | Rasterize attributes of geometries such as "count", "area", or "length" |
| [`as.points`](https://rspatial.github.io/terra/reference/as.points.md)         | Create points from a SpatRaster or SpatVector                           |
| [`as.lines`](https://rspatial.github.io/terra/reference/as.lines.md)           | Create lines from a SpatRaster or SpatVector                            |
| [`as.polygons`](https://rspatial.github.io/terra/reference/as.polygons.md)     | Create polygons from a SpatRaster                                       |
| [`as.contour`](https://rspatial.github.io/terra/reference/contour.md)          | Contour lines from a SpatRaster                                         |
| —————————                                                                      | ——————————————————————————————                                          |

## XXIX. Plotting

### Maps

|                                                                                       |                                                                                         |
|---------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------|
| [`plot`](https://rspatial.github.io/terra/reference/plot.md)                          | Plot a SpatRaster or SpatVector. The main method to create a map                        |
| [`panel`](https://rspatial.github.io/terra/reference/panel.md)                        | Combine multiple plots                                                                  |
| [`points`](https://rspatial.github.io/terra/reference/lines.md)                       | Add points to a map                                                                     |
| [`lines`](https://rspatial.github.io/terra/reference/lines.md)                        | Add lines to a map                                                                      |
| [`polys`](https://rspatial.github.io/terra/reference/lines.md)                        | Add polygons to a map                                                                   |
| [`text`](https://rspatial.github.io/terra/reference/text.md)                          | Add text (such as the values of a SpatRaster or SpatVector) to a map                    |
| [`halo`](https://rspatial.github.io/terra/reference/halo.md)                          | Add text with a halo to a map                                                           |
| [`map.pal`](https://rspatial.github.io/terra/reference/mappal.md)                     | Color palettes for mapping                                                              |
| [`image`](https://rspatial.github.io/terra/reference/image.md)                        | Alternative to plot to make a map with a SpatRaster                                     |
| [`plotRGB`](https://rspatial.github.io/terra/reference/plotRGB.md)                    | Combine three layers (red, green, blue channels) into a single "real color" plot        |
| [`plot<SpatGraticule>`](https://rspatial.github.io/terra/reference/plot_graticule.md) | plot a graticule                                                                        |
| [`sbar`](https://rspatial.github.io/terra/reference/sbar.md)                          | Add a scale bar to a map                                                                |
| [`north`](https://rspatial.github.io/terra/reference/north.md)                        | Add a north arrow to a map                                                              |
| [`inset`](https://rspatial.github.io/terra/reference/inset.md)                        | Add a small inset (overview) map                                                        |
| [`add_legend`](https://rspatial.github.io/terra/reference/legend.md)                  | Add a legend to a map                                                                   |
| [`add_box`](https://rspatial.github.io/terra/reference/box.md)                        | Add a bounding box to a map                                                             |
| [`map_extent`](https://rspatial.github.io/terra/reference/map_extent.md)              | Get the coordinates of a map's axes positions                                           |
| [`dots`](https://rspatial.github.io/terra/reference/dots.md)                          | Make a dot-density map                                                                  |
| [`cartogram`](https://rspatial.github.io/terra/reference/cartogram.md)                | Make a cartogram                                                                        |
| [`persp`](https://rspatial.github.io/terra/reference/persp.md)                        | Perspective plot of a SpatRaster                                                        |
| [`contour`](https://rspatial.github.io/terra/reference/contour.md)                    | Contour plot or filled-contour plot of a SpatRaster                                     |
| [`colorize`](https://rspatial.github.io/terra/reference/RGB.md)                       | Combine three layers (red, green, blue channels) into a single layer with a color-table |
| —————————                                                                             | ——————————————————————————————                                                          |

### Interacting with a map

|                                                                |                                                                           |
|----------------------------------------------------------------|---------------------------------------------------------------------------|
| [`zoom`](https://rspatial.github.io/terra/reference/zoom.md)   | Zoom in to a part of a map by drawing a bounding box on it                |
| [`click`](https://rspatial.github.io/terra/reference/click.md) | Query values of SpatRaster or SpatVector by clicking on a map             |
| [`sel`](https://rspatial.github.io/terra/reference/select.md)  | Select a spatial subset of a SpatRaster or SpatVector by drawing on a map |
| [`draw`](https://rspatial.github.io/terra/reference/draw.md)   | Create a SpatExtent or SpatVector by drawing on a map                     |
| —————————                                                      | ——————————————————————————————                                            |

### Other plots

|                                                                    |                                                                                      |
|--------------------------------------------------------------------|--------------------------------------------------------------------------------------|
| [`plot`](https://rspatial.github.io/terra/reference/plot.md)       | x-y scatter plot of the values of (a sample of) the layers of two SpatRaster objects |
| [`hist`](https://rspatial.github.io/terra/reference/hist.md)       | Histogram of SpatRaster values                                                       |
| [`barplot`](https://rspatial.github.io/terra/reference/barplot.md) | Bar plot of a SpatRaster                                                             |
| [`density`](https://rspatial.github.io/terra/reference/density.md) | Density plot of SpatRaster values                                                    |
| [`pairs`](https://rspatial.github.io/terra/reference/pairs.md)     | Pairs plot for layers in a SpatRaster                                                |
| [`boxplot`](https://rspatial.github.io/terra/reference/boxplot.md) | Box plot of the values of a SpatRaster                                               |
| —————————                                                          | ——————————————————————————————                                                       |

## **Comparison with the raster package**

———————————————————————————————————————

## XXX. New method names

`terra` has a single class `SpatRaster` for which `raster` has three
(`RasterLayer, RasterStack, RasterBrick`). Likewise there is a single
class for vector data `SpatVector` that replaces six `Spatial*` classes.
Most method names are the same, but note the following important
differences in methods names with the `raster` package

|                                  |                                                                                                                                            |
|----------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------|
| **raster package**               | **terra package**                                                                                                                          |
| `raster, brick, stack`           | [`rast`](https://rspatial.github.io/terra/reference/rast.md)                                                                               |
| `rasterFromXYZ`                  | [`rast`](https://rspatial.github.io/terra/reference/rast.md)`( , type="xyz")`                                                              |
| `stack, addLayer`                | [`c`](https://rspatial.github.io/terra/reference/c.md)                                                                                     |
| `addLayer`                       | `add<-`                                                                                                                                    |
| `area`                           | [`cellSize`](https://rspatial.github.io/terra/reference/cellSize.md) or [`expanse`](https://rspatial.github.io/terra/reference/expanse.md) |
| `approxNA`                       | [`approximate`](https://rspatial.github.io/terra/reference/approximate.md)                                                                 |
| `calc`                           | [`app`](https://rspatial.github.io/terra/reference/app.md)                                                                                 |
| `cellFromLine, cellFromPolygon,` | [`cells`](https://rspatial.github.io/terra/reference/cells.md)                                                                             |
| `cellsFromExtent`                | [`cells`](https://rspatial.github.io/terra/reference/cells.md)                                                                             |
| `cellStats`                      | [`global`](https://rspatial.github.io/terra/reference/global.md)                                                                           |
| `clump`                          | [`patches`](https://rspatial.github.io/terra/reference/patches.md)                                                                         |
| `compareRaster`                  | [`compareGeom`](https://rspatial.github.io/terra/reference/compareGeom.md)                                                                 |
| `corLocal`                       | [`focalPairs`](https://rspatial.github.io/terra/reference/focalPairs.md)                                                                   |
| `coordinates`                    | [`crds`](https://rspatial.github.io/terra/reference/crds.md)                                                                               |
| `couldBeLonLat`                  | [`is.lonlat`](https://rspatial.github.io/terra/reference/is.lonlat.md)                                                                     |
| `disaggregate`                   | [`disagg`](https://rspatial.github.io/terra/reference/disaggregate.md)                                                                     |
| `distanceFromPoints`             | [`distance`](https://rspatial.github.io/terra/reference/distance.md)                                                                       |
| `drawExtent, drawPoly, drawLine` | [`draw`](https://rspatial.github.io/terra/reference/draw.md)                                                                               |
| `dropLayer`                      | [`subset`](https://rspatial.github.io/terra/reference/subset.md)                                                                           |
| `extent`                         | [`ext`](https://rspatial.github.io/terra/reference/ext.md)                                                                                 |
| `getValues`                      | [`values`](https://rspatial.github.io/terra/reference/values.md)                                                                           |
| `isLonLat, isGlobalLonLat`       | [`is.lonlat`](https://rspatial.github.io/terra/reference/is.lonlat.md)                                                                     |
| `layerize`                       | [`segregate`](https://rspatial.github.io/terra/reference/segregate.md)                                                                     |
| `layerStats`                     | [`layerCor`](https://rspatial.github.io/terra/reference/layerCor.md)                                                                       |
| `movingFun`                      | [`roll`](https://rspatial.github.io/terra/reference/roll.md)                                                                               |
| `NAvalue`                        | [`NAflag`](https://rspatial.github.io/terra/reference/NAflag.md)                                                                           |
| `nlayers`                        | [`nlyr`](https://rspatial.github.io/terra/reference/dimensions.md)                                                                         |
| `overlay`                        | [`lapp`](https://rspatial.github.io/terra/reference/lapp.md)                                                                               |
| `unstack`                        | [`as.list`](https://rspatial.github.io/terra/reference/as.list.md)                                                                         |
| `projectRaster`                  | [`project`](https://rspatial.github.io/terra/reference/project.md)                                                                         |
| `rasterToPoints`                 | [`as.points`](https://rspatial.github.io/terra/reference/as.points.md)                                                                     |
| `rasterToPolygons`               | [`as.polygons`](https://rspatial.github.io/terra/reference/as.polygons.md)                                                                 |
| `readAll`                        | [`toMemory`](https://rspatial.github.io/terra/reference/toMemory.md)                                                                       |
| `reclassify, subs, cut`          | [`classify`](https://rspatial.github.io/terra/reference/classify.md)                                                                       |
| `sampleRandom, sampleRegular`    | [`spatSample`](https://rspatial.github.io/terra/reference/sample.md)                                                                       |
| `shapefile`                      | [`vect`](https://rspatial.github.io/terra/reference/vect.md)                                                                               |
| `stackApply`                     | [`tapp`](https://rspatial.github.io/terra/reference/tapp.md)                                                                               |
| `stackSelect`                    | [`selectRange`](https://rspatial.github.io/terra/reference/selectRange.md)                                                                 |

## XXXI. Changed behavior

Also note that even if function names are the same in `terra` and
`raster`, their output can be different. In most cases this was done to
get more consistency in the returned values (and thus fewer errors in
the downstream code that uses them). In other cases it simply seemed
better. Here are some examples:

|                                                                                       |                                                                                                                                           |
|---------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|
| [`resample`](https://rspatial.github.io/terra/reference/resample.md)                  | Results are not numerically identical when using `method="bilinear"`, especially at edges, and when going from a high to a low resolution |
| [`as.polygons`](https://rspatial.github.io/terra/reference/as.polygons.md)            | By default, `terra` returns dissolved polygons                                                                                            |
| [`quantile`](https://rspatial.github.io/terra/reference/quantile.md)                  | computes by cell, across layers instead of the other way around                                                                           |
| [`extract`](https://rspatial.github.io/terra/reference/extract.md)                    | By default, `terra` returns a matrix, with the first column the sequential ID of the vectors.                                             |
|                                                                                       | `raster` returns a list (for lines or polygons) or a matrix (for points, but without the ID                                               |
|                                                                                       | column. You can use `list=TRUE` to get the results as a list                                                                              |
| [`values`](https://rspatial.github.io/terra/reference/values.md)                      | `terra` always returns a matrix. `raster` returns a vector for a `RasterLayer`                                                            |
| [`Summary-methods`](https://rspatial.github.io/terra/reference/summarize-generics.md) | With `raster`, `mean(x, y)` and `mean(stack(x, y)` return the same result, a single                                                       |
|                                                                                       | layer with the mean of all cell values. This is also what `terra` returns with                                                            |
|                                                                                       | `mean(c(x, y))`, but with `mean(x, y)` the parallel mean is returned – that is, the                                                       |
|                                                                                       | computation is done layer-wise, and the number of layers in the output is the same as                                                     |
|                                                                                       | that of `x` and `y` (or the larger of the two if they are not the same). This affects                                                     |
|                                                                                       | all summary functions (`sum`, `mean`, `median`, `which.min`, `which.max`, `min`, `max`,                                                   |
|                                                                                       | `prod`, `any`, `all`, `stdev`), except `range`, which is not implemented for this case                                                    |
|                                                                                       | (you can use `min` and `max` instead)                                                                                                     |
| —————————                                                                             | ——————————————————————————————                                                                                                            |

## Contributors

Except where indicated otherwise, the methods and functions in this
package were written by Robert Hijmans. The configuration scripts were
written by Roger Bivand. Some of code using the GEOS library was adapted
from code by Edzer Pebesma for `sf`. Emanuele Cordano contributed
functionality for catchment related computations. Andrew Gene Brown,
Márcia Barbosa, Michael Chirico, Krzysztof Dyba, Barry Rowlingson, and
Michael D. Sumner also made important contributions

This package is an attempt to climb on the shoulders of giants (GDAL,
PROJ, GEOS, NCDF, GeographicLib, Rcpp, R). Many people have contributed
by asking questions or [raising
issues](https://github.com/rspatial/terra). Feedback and suggestions by
Kendon Bell, Jean-Luc Dupouey, Sarah Endicott, Derek Friend, Alex Ilich,
Agustin Lobo, Gerald Nelson, Jakub Nowosad, and Monika Tomaszewska have
been especially helpful.

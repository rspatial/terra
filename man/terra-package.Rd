\name{terra-package}
\alias{terra-package}

\alias{terra}

\docType{package}

\title{The terra package}

\description{

\code{terra} provides methods to manipulate geographic (spatial) data in "raster" and "vector" form. Raster data divides space into rectangular cells (pixels) and it is commonly used to represent spatially continuous phenomena, such as elevation or the weather. Satellite images also have this data structure. In contrast, "vector" spatial data (points, lines, polygons) is typically used to represent discrete spatial entities. 

The package implements two classes (R data types) to reprsent these data types: \code{SpatRaster} and \code{SpatVector}. \code{SpatRaster} supports handling large raster files that cannot be loaded into memory; local, focal, zonal, and global raster operations; polygon, line and point to raster conversion; integration with modeling methods to make spatial predictions; and more. \code{SpatVector} supports all types of geometric operations such as intersections.

Additional classes include \code{SpatExtent}, which is used to define a spatial extent (bounding box); \code{SpatRasterDataset} to represents a collection of sub-datasets for the same area. Each subdataset is a SpatRaster with possibly many layers, and may, for example, represent different weather variables; and \code{SpatRasterCollection} and \code{SpatVectorCollection} that are vectors of SpatRaster or SpatVector.

The \code{terra} package is conceived as a replacement of the \code{raster} package. \code{terra} has a very similar, but simpler, interface, and it is faster than \code{raster}. At the bottom of this page there is a table that shows differences in the methods between the two packages. 

Below is a list of some of the most important methods grouped by theme. Some of these may not have been implemented yet (they are not hyperlinked). 
}

\section{-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}{

\bold{SpatRaster}

---------------------------------------------------------------------------------------------------------------------}

\section{I. Creating, combining and subsetting SpatRaster objects}{

  \tabular{ll}{
    \code{\link{rast}}\tab Create a SpatRaster from scratch, file, or another object\cr
    \code{\link{c}} \tab Combine SpatRasters (multiple layers)\cr
    \code{\link{add<-}} \tab Add a SpatRaster to another one\cr
    \code{\link{subset}} or \code{[[}, or \code{$} \tab Select layers of a SpatRaster\cr
    \code{\link{selectRange}} \tab Select cell values from different layers using an index layer\cr
--------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}


\section{II. Changing the spatial extent and/or resolution of a SpatRaster}{
Also see the methods in section VIII

  \tabular{ll}{
    \code{\link{merge}} \tab Combine SpatRasters with different extents (but same origin and resolution) \cr
    \code{\link{mosaic}} \tab Combine SpatRasters with different extents and a function for the values in overlapping areas \cr
    \code{\link{crop}} \tab Select a geographic subset of a SpatRaster \cr
    \code{\link{expand}} \tab Enlarge a SpatRaster \cr
    \code{\link{trim}} \tab Trim a SpatRaster by removing exterior rows and/or columns that only have NAs\cr
    \code{\link{aggregate}} \tab Combine cells of a SpatRaster to create larger cells \cr
    \code{\link{disaggregate}} \tab Subdivide cells \cr
    \code{\link{resample}} \tab Resample (warp) values to a SpatRaster with a different origin and/or resolution \cr 
	\code{\link{project}} \tab Project (warp) values to a SpatRaster with a different coordinate reference system \cr
    \code{\link{shift}} \tab Adjust the location of SpatRaster \cr
    \code{\link{flip}} \tab  Flip values horizontally or vertically \cr
    \code{\link{rotate}} \tab Rotate values around the date-line (for lon/lat data) \cr
    \code{\link{t}} \tab Transpose a SpatRaster\cr
--------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}

  
\section{III. Local (cell based) computation}{

  \tabular{ll}{
    \code{\link{app}} \tab Apply a function to cells, across layers (as in \code{base::apply}) \cr
    \code{\link{tapp}} \tab Apply a function to groups of layers (as in \code{base::tapply})\cr
    \code{\link{lapp}} \tab Apply a function using the layers as variables\cr
    \code{\link{rapp}} \tab Apply a function to a spatially variable range of layers\cr
    \code{\link{arith}} \tab Use a function to compute new values for all cells \cr
    \code{\link{Arith-methods}} \tab Standard arith methods (\code{+, -, *, ^, \%\%, \%/\%, /}) \cr
	\code{\link{Math-methods}} \tab Math methods like \code{abs, sqrt, trunc, log, log10, exp, sin, round} \cr	
    \code{\link{Logic-methods}} \tab Boolean methods (\code{!, &, |}) \cr
    \code{\link{Summary-methods}} \tab Summary methods (\code{mean, max, min, median, sum, range, prod, any, all, stdev, which.min, which.max}) \cr
    \code{\link{Compare-methods}} \tab Comparison methods (\code{==, !=, >, <, <=, >=}) \cr
    \code{\link{area}} \tab Compute the area of cells \cr
    \code{\link{classify}} \tab (Re-)classify values \cr
    \code{\link{cover}} \tab First layer covers second layer except where the first layer is \code{NA} \cr
    \code{\link{init}} \tab Initialize cells with new values \cr	
    \code{\link{mask}} \tab Use values from first SpatRaster except where cells of the mask SpatRaster are \code{NA} (or another value of choice)\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}

\section{IV. Zonal and global computation}{
  \tabular{ll}{
    \code{\link{area}} \tab Compute the total area covered by cells \cr
    \code{\link{crosstab}} \tab Cross-tabulate two SpatRasters\cr
    \code{\link{freq}} \tab Frequency table of SpatRaster cell values \cr
    \code{\link{global}} \tab Summarize SpatRaster cell values with a function \cr
    \code{\link{quantile}} \tab Quantiles \cr
    \code{\link{stretch}} \tab Stretch values \cr
    \code{\link{scale}} \tab Scale values \cr
    \code{\link[terra]{summary}} \tab Summary of the values of a SpatRaster (quartiles and mean) \cr
    \code{\link{unique}} \tab Get the unique values in a SpatRaster \cr
    \code{\link{zonal}} \tab Summarize a SpatRaster by zones in another SpatRaster \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}


\section{V. Focal and other spatial contextual computation}{

  \tabular{ll}{
    \code{\link{focal}} \tab Focal (neighborhood; moving window) functions \cr
    \code{\link{adjacent}} \tab Identify cells that are adjacent to a set of cells of a SpatRaster \cr
    \code{\link{boundaries}} \tab Detection of boundaries (edges)\cr
    \code{\link{distance}} \tab Shortest distance to a cell that is not \code{NA} or to or from a vector object\cr
    \code{direction} \tab Direction (azimuth) to or from cells that are not \code{NA}\cr
    \code{localFun} \tab Local association (using neighborhoods) functions \cr
    \code{patches} \tab Find patches \cr
    \code{\link{terrain}} \tab Compute slope, aspect and other terrain characteristics from elevation data \cr	
    \code{\link{autocor}} \tab Compute global or local spatial autocorrelation  \cr	
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}

\section{VI. Model predictions}{

  \tabular{ll}{
    \code{\link{predict}} \tab Predict a non-spatial model to a SpatRaster \cr
    \code{\link{interpolate}} \tab Predict a spatial model to a SpatRaster \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}


\section{VII. Accessing cell values}{

Apart from the function listed below, you can also use indexing with \code{[} with cell numbers, and row and/or column numbers \cr

\tabular{ll}{
    \code{\link{values}} \tab Get or set all cell values (fails with very large rasters)\cr
    \code{\link{setValues}} \tab Set new values to the cells of a SpatRaster \cr
    \code{\link{as.matrix}} \tab Get cell values as a matrix \cr
    \code{\link{as.array}} \tab Get cell values as an array \cr
    \code{\link{extract}} \tab Extract cell values from a SpatRaster (e.g., by cell, coordinates, polygon)\cr
    \code{\link{spatSample}} \tab Regular or random sample \cr
    \code{\link{minmax}} \tab Get the minimum and maximum value of the cells of a SpatRaster (if known) \cr
    \code{\link{setMinMax}} \tab Compute the minimum and maximum value of a SpatRaster if these are not known \cr
    \code{\link{extract}} \tab spatial queries of a SpatRaster with a SpatVector\cr
	
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}


\section{VIII. Getting and setting SpatRaster dimensions }{

Get or set basic parameters of SpatRasters. If there are values associated with a SpatRaster object (either in memory or via a link to a file) these are lost when you change the number of columns or rows or the resolution. This is not the case when the extent is changed (as the number of columns and rows will not be affected). Similarly, with \bold{crs} you can set the coordinate reference system, but this does not transform the data (see \link{project} for that).   

  \tabular{ll}{
    \code{\link{ncol}}\tab The number of columns \cr
    \code{\link{nrow}} \tab The number of rows \cr
    \code{\link{ncell}} \tab The number of cells (can not be set directly, only via ncol or nrow) \cr
    \code{\link{res}} \tab The resolution (x and y) \cr
    \code{\link{nlyr}} \tab Get or set the number of layers \cr	
    \code{\link{names}} \tab Get or set the layer names \cr	
    \code{\link{xres}} \tab The x resolution (can be set with res) \cr
    \code{\link{yres}} \tab The y resolution (can be set with res)\cr
    \code{\link{xmin}} \tab The minimum x coordinate (or longitude) \cr
    \code{\link{xmax}} \tab The maximum x coordinate (or longitude) \cr
    \code{\link{ymin}} \tab The minimum y coordinate (or latitude) \cr
    \code{\link{ymax}} \tab The maximum y coordinate (or latitude) \cr
	\code{\link{ext}} \tab Get or set the extent (minimum and maximum x and y coordinates (a.k.a. "bounding box") \cr
	\code{\link{origin}} \tab The origin of a SpatRaster\cr
    \code{\link{crs}} \tab The coordinate reference system (map projection) \cr	
    \code{\link{is.lonlat}} \tab Test if an object has (or may have) a longitude/latitude coordinate reference system; and if it has global coverage \cr
    \code{\link{sources}} \tab Get the filename(s) to which a SpatRaster is linked \cr	
    \code{\link{compareGeom}} \tab Compare the geometry of SpatRasters \cr		
    \code{\link{NAflag}} \tab Set the \code{NA} value (for reading from a file with insufficient metadata) \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}


\section{IX. Computing row, column, cell numbers and coordinates}{

Cell numbers start at 1 in the upper-left corner. They increase within rows, from left to right, and then row by row from top to bottom. Likewise, row numbers start at 1 at the top of the raster, and column numbers start at 1 at the left side of the raster. 

  \tabular{ll}{
    \code{\link{xFromCol}} \tab x-coordinates from column numbers \cr
    \code{\link{yFromRow}} \tab y-coordinates from row numbers \cr
    \code{\link{xFromCell}} \tab x-coordinates from row numbers \cr
    \code{\link{yFromCell}} \tab y-coordinates from cell numbers \cr
    \code{\link{xyFromCell}} \tab x and y coordinates from cell numbers \cr
    \code{\link{colFromX}} \tab Column numbers from x-coordinates (or longitude) \cr
    \code{\link{rowFromY}} \tab Row numbers from y-coordinates (or latitude) \cr
    \code{\link{rowColFromCell}} \tab Row and column numbers from cell numbers\cr
    \code{\link{cellFromXY}} \tab Cell numbers from x and y coordinates \cr
    \code{\link{cellFromRowCol}} \tab Cell numbers from row and column numbers \cr
    \code{\link{cellFromRowColCombine}} \tab Cell numbers from all combinations of row and column numbers \cr    
	\code{\link{cells}}\tab Cell numbers from an SpatVector or SpatExtent\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}


\section{X. Writing SpatRaster files}{

  \tabular{ll}{
    Basic\tab\cr
    \code{\link{writeRaster}} \tab Write all values of SpatRaster to disk \cr
   .\cr
    \code{\link{writeCDF}} \tab Write netCDF files \cr
   .\cr

   Advanced\tab\cr
   \code{\link{blockSize}} \tab Get suggested block size for reading and writing \cr
   \code{\link{writeStart}} \tab Open a file for writing \cr
   \code{\link{writeValues}} \tab Write some values  \cr
   \code{\link{writeStop}} \tab Close the file after writing \cr
    --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}	

\section{XI. Miscellaneous SpatRaster methods}{

\tabular{ll}{
    \code{\link{terraOptions}} \tab Show, set, save or get session options \cr
    \code{\link{sources}} \tab Show the data sources of a SpatRaster \cr
    \code{\link{tmpFiles}} \tab Show or remove temporary files \cr
    \code{\link{canProcessInMemory}} \tab Test whether a file can be created in memory \cr
    \code{\link{readStart}} \tab Open file connections for efficient multi-chunck reading \cr
    \code{\link{readStop}} \tab Close file connections \cr
    \code{\link{inMemory}} \tab Are the cell values in memory? \cr
    \code{\link{fromDisk}} \tab Are the cell values read from a file? \cr	
}
}


\section{-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}{

\bold{SpatRasterDataSet}

---------------------------------------------------------------------------------------------------------------------}


\section{XII. SpatRasterDataset}{

A SpatRasterDataset contains SpatRaster objects that are sub-datasets for the same area. They all have the same extent and resolution.

\tabular{ll}{
  \code{\link{sds}} \tab Create a SpatRasterDataset \cr
  \code{[} or \code{$} \tab Extract a SpatRaster \cr
  \code{\link{names}} \tab Get the names of the sub-datasets \cr
}
}


\section{-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}{

\bold{SpatVector}

---------------------------------------------------------------------------------------------------------------------}


\section{XIII. Create and combine SpatVector objects}{
\tabular{ll}{
  \code{\link{vect}} \tab Create a SpatRaster from a file (e.g. shapefile) or from another object\cr
  \code{\link{c}} \tab \bold{append} vector data of the same (vector) type ("rbind") \cr
  \code{\link{unique}} \tab remove duplicates \cr
  \code{\link{project}} \tab Project a SpatVector to a different coordinate reference system \cr
  \code{\link{writeVector}} \tab Write SpatVector data to disk \cr
  \code{\link{centroids}} \tab Get the centroids of a SpatVector\cr
  \code{\link{voronoi}} \tab Voronoi diagram \cr
  \code{\link{delauny}} \tab Delauny triangles\cr
  \code{\link{convexhull}} \tab Compute the convex hull of a SpatVector \cr  
  \code{\link{fill}}\tab Remove or extract holes from polygons\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr  
}}


\section{XIV. Properties of SpatVector objects}{
\tabular{ll}{
  \code{\link{geom}} \tab returns the geometries as matrix or WKT\cr
  \code{\link{linearUnits}} \tab returns the linear units of the crs (in meter)\cr
  \code{\link{ncol}}\tab The number of columns (of the attributes)\cr
  \code{\link{nrow}} \tab The number of rows (of the geometries and attributes)\cr
  \code{\link{names}} \tab Get or set the layer names \cr	
  \code{\link{ext}} \tab Get the extent (minimum and maximum x and y coordinates (a.k.a. "bounding box") \cr
  \code{\link{crs}} \tab The coordinate reference system (map projection) \cr	
  \code{\link{is.lonlat}} \tab Test if an object has (or may have) a longitude/latitude coordinate reference system; and if it has global coverage \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr  
}}


\section{XV. Geometric queries}{
\tabular{ll}{
  \code{\link{relate}} \tab geometric relationships such as "intersects", "overlaps", and "touches"\cr
  \code{\link{adjacent}} \tab find adjacent polygons\cr
  \code{\link{near}} \tab find nearby geometries\cr
  \code{\link{area}} \tab computes the area covered by polygons\cr
  \code{\link{perimeter}} \tab computes the length of the perimeter of polygons, and the length of lines\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr  
}}



\section{XVI. Geometric operations}{

The name in \bold{bold} is the equivalent command in ArcGIS. 

\tabular{ll}{
  \code{\link{erase}} or "-" \tab \bold{erase} (parts of) geometries\cr
  \code{\link{intersect}} or "*"  \tab \bold{intersect} geometries\cr
  \code{\link{union}} or "+" \tab \bold{union}. Merge geometries\cr
  \code{\link{cover}} \tab \bold{update} and \bold{identity} a polygons\cr
  \code{\link{symdif}} \tab\bold{symmetrical difference} of two polygons \cr 
  \code{\link{aggregate}} \tab \bold{dissolve} smaller polygons into larger ones \cr
  \code{\link{buffer}} \tab \bold{buffer} \cr
  \code{\link{disaggregate}} \tab \bold{explode}: turn multi-geometries into separate geometries \cr
  \code{\link{crop}} \tab \bold{clip} geometries using a rectangle (SpatExtent)\cr
   .\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr  
}}


\section{XVII. SpatVector attributes}{

\tabular{ll}{
  \code{\link{extract}} \tab spatial queries between SpatVector and SpatVector (e.g. point in polygons) \cr
  \code{\link{select}} \tab \bold{select} - interactively select geometries\cr
  \code{\link{click}} \tab \bold{identify} attributes by clicking on a map\cr
  \code{\link{merge}} \tab \bold{Join table} \cr
  \code{\link{as.data.frame}} \tab get attributes as a data.frame\cr
}}


\section{XVIII. Change geometries (for display, experimentation)}{

\tabular{ll}{
  \code{\link{shift}} \tab change the position geometries \cr
  \code{\link{spin}}\tab rotate geometries around an origin\cr
  \code{\link{rescale}} \tab shrink (or expand) geometries, for example to make an inset map \cr
  \code{\link{flip}} \tab flip geometries vertially or horizontally\cr
  \code{\link{t}} \tab transpose geometries (switch x and y)\cr
}}



\section{-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}{

\bold{Spat* Collections}

---------------------------------------------------------------------------------------------------------------------}


\section{XIX. Collections}{
A SpatRasterCollection is a vector of SpatRaster objects. Unlike for a SpatRasterDataset, there the extent and resolution of the elements do not need to match each other. A SpatVectorCollection is a vector of SpatVector objects. 

\tabular{ll}{
\code{\link{svc}} \tab create a SpatRasterCollection\cr
\code{\link{length}} \tab how many elements does the collection have?\cr
\code{\link{[}} \tab extract an element\cr
}}



\section{-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}{

\bold{SpatExtent}

---------------------------------------------------------------------------------------------------------------------}



\section{XX. SpatExtent}{

  \tabular{ll}{
    \code{\link{extent}} \tab Create a SpatExtent object \cr
    \code{\link{intersect}} \tab Intersect two SpatExtent objects, same as \code{-}  \cr
    \code{\link{union}} \tab Combine two SpatExtent objects, same as \code{+} \cr
    \code{\link{Math-methods}} \tab round/floor/ceiling of a SpatExtent \cr			
    \code{\link{align}} \tab Align a SpatExtent with a SpatRaster \cr
    \code{\link{draw}} \tab Create a SpatExtent by drawing it on top of a map (plot) \cr		
  }
}




\section{-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}{

\bold{General methods}

---------------------------------------------------------------------------------------------------------------------}


\section{XXI. Data type conversion}{

You can coerce SpatRasters to Raster* objects after loading the \code{raster} package with \code{as(object, "Raster")}, or \code{raster(object)} or \code{brick(object)} or \code{stack(object)}
  
  \tabular{ll}{
    \code{\link{rast}} \tab SpatRaster from matrix and other objects\cr
    \code{\link{rasterize}} \tab Rasterizing points, lines or polygons\cr
    \code{\link{as.points}} \tab Create points from a SpatRaster or SpatVector \cr
    \code{\link{as.lines}} \tab Create points from a SpatRaster or SpatVector\cr
    \code{\link{as.polygons}} \tab Create polygons from a SpatRaster \cr
    \code{\link{as.contour}}  \tab Contour lines from a SpatRaster \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}

\section{XXII. Plotting}{

  \tabular{ll}{
    Maps \tab \cr
    \code{\link{plot}} \tab Plot a SpatRaster or SpatVector. The main method to create a map \cr
	\code{\link{points}} \tab Add points to a map\cr
	\code{\link{lines}} \tab Add lines to a map\cr
	\code{\link{polys}} \tab Add polygons to a map \cr
    \code{\link{plotRGB}} \tab Combine three layers (red, green, blue channels) into a single "real color" image \cr
    \code{\link{image}} \tab Alternative way to plot a SpatRaster \cr
    \code{\link{persp}} \tab Perspective plot of a SpatRaster \cr
    \code{\link{contour}} \tab Contour plot or filled-contour plot of a SpatRaster \cr
    \code{\link{text}} \tab Plot the values of a SpatRaster or SpatVector on top of a map \cr
    \code{\link{inset}} \tab Add a small inset (overview) map\cr
    \code{\link{sbar}} \tab Add a scalebar\cr
    .\cr
    Interacting with a map \tab \cr
    \code{\link{zoom}} \tab Zoom in to a part of a map \cr
    \code{\link{click}} \tab Query values of SpatRaster or SpatVector by clicking on a map \cr
    \code{\link{select}} \tab Select a spatial subset of a SpatRaster or SpatVector \cr
    \code{\link{draw}} \tab Create a SpatExtent or SpatVector by drawing on a map \cr
    .\cr
	
	Other plots\tab \cr
    \code{\link{plot}} \tab x-y scatter plot of the values of two SpatRaster objects\cr
    \code{\link{hist}} \tab Histogram of SpatRaster values \cr
    \code{\link{barplot}} \tab Barplot of a SpatRaster \cr
    \code{\link{density}} \tab Density plot of SpatRaster values \cr
    \code{\link{pairs}} \tab Pairs plot for layers in a SpatRaster \cr
    \code{\link{boxplot}} \tab Box plot of the values of a SpatRaster\cr
  }
}


\section{XXIII. Comparison with the raster package}{

\code{terra} has a single class \code{SpatRaster} for which \code{raster} has three (\code{RasterLayer, RasterStack, RasterBrick}). Likewise there is a single class for vector data \code{SpatVector} that replaces six \code{Spatial*} classes. Most method names are the same, but 
note the following important differences in methods names with the 'raster' package

\tabular{ll}{
\bold{raster package} \tab \bold{terra package}\cr
\code{raster, brick, stack}\tab\code{\link{rast}}\cr
\code{rasterFromXYZ}\tab\code{\link{rast}( , type="xyz")}\cr
\code{stack, addLayer} (combining Raster* objects or files)\tab\code{\link{c}}\cr
\code{addLayer} \tab\code{\link{add<-}}\cr
\code{extent}\tab\code{\link{ext}}\cr
\code{calc}\tab\code{\link{app}} and \code{\link{arith}}\cr
\code{overlay}\tab\code{\link{lapp}}\cr
\code{stackApply}\tab\code{\link{tapp}}\cr
\code{extend}\tab\code{\link{expand}}\cr
\code{nlayers}\tab\code{\link{nlyr}}\cr
\code{NAvalue}\tab\code{\link{NAflag}}\cr
\code{stackSelect}\tab\code{\link{selectRange}}\cr
\code{reclassify, subs, cut}\tab\code{\link{classify}}\cr
\code{cellStats}\tab\code{\link{global}}\cr
\code{projectRaster}\tab\code{\link{project}}\cr
\code{dropLayer}\tab\code{\link{subset}}\cr
\code{isLonLat, isGlobalLonLat, couldBeLonLat}\tab\code{\link{is.lonlat}}\cr
\code{shapefile}\tab\code{\link{vect}}\cr
\code{gridDistance, distanceFromPoints} \tab \code{\link{distance}}\cr
\code{drawExtent, drawPoly, drawLine} \tab \code{\link{draw}}\cr
\code{compareRaster}\tab\code{\link{compareGeom}}\cr
\code{sampleRandom, sampleRegular}\tab\code{\link{spatSample}}\cr
\code{rasterToPoints}\tab\code{\link{as.points}}\cr
\code{rasterToPolygons}\tab\code{\link{as.polygons}}\cr
\code{cellFromLine, cellFromPolygon, cellsFromExtent}\tab\code{\link{cells}}\cr
\code{layerize}\tab\code{\link{separate}}\cr
\code{clump}\tab\code{patches}\cr
}

Also note that even if function names are the same in \code{terra} and \code{raster}, their output can be different. In most cases to get more consistency in the returned values (and thus fewer errors in the downstream code that uses them). It other cases it simply seemed better. Here are some examples:

\tabular{ll}{
\code{\link{area}}\tab By default, \code{terra} returns the summed area of the raster cells that are not \code{NA}. \code{raster} returns a \code{RasterLayer}. Use \code{sum=FALSE} in \code{terra} for that behavior\cr
- \tab \cr
\code{\link{as.polygons}}\tab By default, \code{terra} returns dissolved polygons\cr
- \tab \cr
\code{\link{extract}}\tab By default, \code{terra} returns a matrix, with the first column the sequential ID of the vectors. \code{raster} returns a list (for lines or polygons) or a matrix (for points, but without the ID column. You can use \code{list=TRUE} to get the results as a list\cr
- \tab \cr
\code{\link{values}}\tab \code{terra} always returns a matrix. \code{raster} returns a vector for a \code{RasterLayer}\cr
- \tab \cr
\code{\link{Summary-methods}}\tab With \code{raster}, \code{mean(x, y)} and \code{mean(stack(x, y)} return the same result, a single layer with the mean of all cell values. That is what \code{terra} returns with \code{mean(c(x, y))}, but with \code{mean(x, y)} the parallel mean is returned -- that is, the computation is done layer-wise, and the number of layers in the output is the same as that of \code{x} and \code{y} (or the larger of the two if they are not the same). This affects all summary functions (\code{sum}, \code{mean}, \code{median}, \code{which.min}, \code{which.max}, \code{min}, \code{max}, \code{prod}, \code{any}, \code{all}, \code{stdev}), except \code{range}, which is not implemented for this case (you can use \code{min} and \code{max} instead) \cr
}
}


\section{-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}{.}


\author{Except where indicated otherwise, the methods and functions in this package were written by Robert Hijmans. The configuration scripts were written by Roger Bivand. Some of the underlying C++ code for GDAL/GEOS was adapted from code written by Edzer Pebesma for \code{sf}. The progress bar code is by Karl Forner (RcppProgress). Jeroen Ooms provided the compiled GDAL and GEOS libraries for installation on windows}

\section{Acknowledgments}{
This package is an attempt to climb on the shoulders of giants (GDAL, PROJ, GEOS, NCDF, GeographicLib, Rcpp, R). Many people have contributed by asking questions or \href{https://github.com/rspatial/terra}{filing bug reports}. The feedback and suggestsions by Kendon Bell, Jean-Luc Dupouey, Gerald Nelson, and Michael Sumner have been especially helpful.
}

\keyword{ package }
\keyword{ spatial }


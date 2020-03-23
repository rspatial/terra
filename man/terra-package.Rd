\name{terra-package}
\alias{terra-package}

\alias{terra}

\docType{package}

\title{The "terra" package}

\description{
The "terra" package implements classes for spatial data (see \link{SpatRaster-class}) and supports handling large raster files that cannot be loaded into memory; local, focal, zonal, and global raster operations; polygon, line and point to raster conversion; integration with modeling methods to make spatial predictions; and more.

The package is conceived as a replacement of the "raster" package. It has a very similar, but simpler, interface, and it is faster. Like the "raster" package, "terra" provides classes and functions to manipulate geographic (spatial) data in "raster" and "vector" data types. Raster data divides space into cells (rectangles; pixels) of equal size (in units of the coordinate reference system). Such continuous spatial data are also referred to as "grid" data, and be contrasted with discrete "vector" spatial data such as points, lines, polygons. A major simplification is that "terra" has a single class for raster data, "SpatRaster", whereas the "raster" package has three (RasterLayer, RasterStack, RasterBrick). Likewise there is a single class for vector data ("SpatVector"). Like "raster" this package should be particularly useful when using very large datasets that can not be loaded into the computer's memory. Functions will work correctly, because they they process large files in chunks, i.e., they read, compute, and write blocks of data, without loading all values into memory at once. 

Note the following important differences in function names with the raster package

\tabular{ll}{
raster package \tab terra package\cr
\code{raster, brick}\tab\code{\link{rast}}\cr
\code{stack} (for creating an object from file)\tab\code{\link{rast}}\cr
\code{stack} (for combining Raster* objects)\tab\code{\link{c}}\cr
\code{rasterFromXYZ}\tab\code{\link{rast}( , type="xyz")}\cr
\code{extent}\tab\code{\link{ext}}\cr
\code{calc}\tab\code{\link{app}}\cr
\code{stackApply}\tab\code{\link{tapp}}\cr
\code{reclassify, subs, cut}\tab\code{\link{classify}}\cr
\code{cellStats}\tab\code{\link{global}}\cr
\code{projectRaster}\tab\code{\link{warp}, \link{project}}\cr
\code{shapefile}\tab\code{\link{vect}}\cr
\code{gridDistance, distanceFromPoints} \tab \code{\link{distance}}\cr
\code{drawExtent, drawPoly, drawLine} \tab \code{\link{draw}}\cr
\code{nlayers}\tab\code{\link{nlyr}}\cr
\code{stackSelect}\tab\code{\link{collapse}}\cr
\code{compareRaster}\tab\code{\link{compareGeom}}\cr
\code{sampleRegular}\tab\code{\link{spatSample}}\cr
}

Below is a list of some of the most important methods grouped by theme. Some of these may not have been implemented yet. 
}


\section{I. Creating SpatRaster objects}{

SpatRaster objects can be created, from scratch, files, or from objects of other classes, with the following functions:
  \tabular{ll}{
    \code{\link{rast}}\tab Create a SpatRaster from scratch, file, or another object\cr
    \code{\link{c}} \tab Combine SpatRasters (multiple layers) (like raster::stack)\cr
    \code{\link{subset}} \tab Select layers of a SpatRaster\cr
--------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}


\section{II. Changing the spatial extent and/or resolution of a SpatRaster}{

  \tabular{ll}{
    \code{\link{merge}} \tab Combine SpatRasters with different extents (but same origin and resolution) \cr
    \code{mosaic} \tab Combine SpatRasters with different extents and a function for the values in overlapping areas \cr
    \code{\link{crop}} \tab Select a geographic subset of a SpatRaster \cr
    \code{\link{extend}} \tab Enlarge a SpatRaster \cr
    \code{\link{trim}} \tab Trim a SpatRaster by removing exterior rows and/or columns that only have NAs\cr
    \code{\link{aggregate}} \tab Combine cells of a SpatRaster to create larger cells \cr
    \code{\link{disaggregate}} \tab Subdivide cells \cr
    \code{\link{warp}} \tab Warp values to a SpatRaster with a different origin and/or resolution and/or coordinate reference system \cr	
    \code{\link{shift}} \tab Adjust the location of SpatRaster \cr
    \code{\link{flip}} \tab  Flip values horizontally or vertically \cr
    \code{\link{rotate}} \tab Rotate values around the date-line (for lon/lat data) \cr
    \code{\link{t}} \tab Transpose a SpatRaster\cr
--------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}

  
\section{III. Local (cell based) computation}{

  \tabular{ll}{
    \code{\link{Arith-methods}} \tab Arith functions (\code{+, -, *, ^, \%\%, \%/\%, /}) \cr
	\code{\link{Math-methods}} \tab Math functions like \code{abs, sqrt, trunc, log, log10, exp, sin, round} \cr	
    \code{\link{Logic-methods}} \tab Logic functions (\code{!, &, |}) \cr
    \code{\link{Summary-methods}} \tab Summary functions (\code{mean, max, min, range, prod, sum, any, all}) \cr
    \code{\link{Compare-methods}} \tab Compare functions (\code{==, !=, >, <, <=, >=}) \cr
    \code{\link{app}} \tab Apply a function to cells of a SpatRaster (as in \code{base::apply}) \cr
    \code{\link{tapp}} \tab Apply a function to cells of a SpatRaster by groups of layers (as in \code{base::tapply})\cr
    \code{\link{overlay}} \tab Computations on multiple SpatRaster objects \cr
    \code{\link{cover}} \tab First layer covers second layer except where the first layer is \code{NA} \cr
    \code{\link{mask}} \tab Use values from first SpatRaster except where cells of the mask SpatRaster are \code{NA} (or another value of choice)\cr
    \code{\link{classify}} \tab (Re-)classify values \cr
    \code{\link{init}} \tab Initialize cells with new values \cr	
    \code{\link{area}} \tab Compute area of cells (for longitude/latitude data) \cr
    \code{\link{collapse}} \tab Select cell values from different layers using an index layer\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}

\section{IV. Zonal and global computation}{
  \tabular{ll}{
    \code{\link{zonal}} \tab Summarize a SpatRaster by zones in another SpatRaster \cr
    \code{\link{global}} \tab Summarize SpatRaster cell values with a function \cr
    \code{\link{unique}} \tab Get the unique values in a SpatRaster \cr
    \code{\link{freq}} \tab Frequency table of SpatRaster cell values \cr
    \code{crosstab} \tab Cross-tabulate two SpatRasters\cr
    \code{quantile} \tab Quantiles \cr
    \code{summary} \tab Summary of the values of a SpatRaster (quartiles and mean) \cr
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
    \code{clump} \tab Find clumps (patches) \cr
    \code{terrain} \tab Compute slope, aspect and other terrain characteristics from elevation data \cr	
    \code{autocor} \tab Compute global or local Moran or Geary indices of spatial autocorrelation  \cr	
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

\section{VII. Data type conversion}{

You can coerce SpatRasters to Raster* objects after loading the \code{raster} package, using \code{as}, as in \code{as(object, "RasterBrick")}, or \code{raster(object)} or \code{brick(object)} or \code{stack(object)}
  
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

\section{VIII. Accessing values of SpatRaster cells}{

Apart from the function listed below, you can also use indexing with \code{[} for cell numbers, and \code{[[} for row / column number combinations \cr

\tabular{ll}{
    \code{\link{values}} \tab Get all cell values (fails with very large rasters), or a row of values (safer) \cr
    \code{\link{as.matrix}} \tab Get cell values as a matrix \cr
    \code{\link{as.array}} \tab Get cell values as an array \cr
    \code{\link{extract}} \tab Extract cell values from a SpatRaster (e.g., by cell, coordinates, polygon)\cr
    \code{sampleRandom} \tab Random sample \cr
    \code{\link{spatSample}} \tab Regular sample \cr
    \code{\link{minmax}} \tab Get the minimum and maximum value of the cells of a SpatRaster (if known) \cr
    \code{\link{setMinMax}} \tab Compute the minimum and maximum value of a SpatRaster if these are not known \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}


\section{IX. Plotting}{

  \tabular{ll}{
    \bold{Maps} \tab \cr
    \code{\link{plot}} \tab Plot a SpatRaster. The main method to create a map \cr
    \code{\link{plotRGB}} \tab Combine three layers (red, green, blue channels) into a single "real color" image \cr
    \code{\link{image}} \tab Plot a SpatRaster with the image function \cr
    \code{\link{persp}} \tab Perspective plot of a SpatRaster \cr
    \code{\link{contour}} \tab Contour plot or filled-contour plot of a SpatRaster \cr
    \code{\link{text}} \tab Plot the values of a SpatRaster or SpatVector on top of a map \cr
    .\cr
    \bold{Interacting with a map} \tab \cr
    \code{\link{zoom}} \tab Zoom in to a part of a map \cr
    \code{\link{click}} \tab Query values of SpatRaster or SpatVector by clicking on a map \cr
    \code{\link{select}} \tab Select a spatial subset of a SpatRaster or SpatVector \cr
    \code{\link{draw}} \tab Create a SpatExtent or SpatVector by drawing on a map \cr
    .\cr
	
	\bold{Other plots} \tab \cr
    \code{\link{plot}} \tab x-y scatter plot of the values of two SpatRaster objects\cr
    \code{\link{hist}} \tab Histogram of SpatRaster values \cr
    \code{\link{barplot}} \tab Barplot of a SpatRaster \cr
    \code{\link{density}} \tab Density plot of SpatRaster values \cr
    \code{\link{pairs}} \tab Pairs plot for layers in a SpatRaster \cr
    \code{\link{boxplot}} \tab Box plot of the values of a SpatRaster\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}

\section{X. Getting and setting SpatRaster dimensions }{

Basic parameters of existing SpatRasters can be obtained, and in most cases changed. If there are values associated with a SpatRaster object (either in memory or via a link to a file) these are lost when you change the number of columns or rows or the resolution. This is not the case when the extent is changed (as the number of columns and rows will not be affected). Similarly, with \bold{projection} you can set the projection, but this does not transform the data (see \link{warp} for that).   

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
    \code{\link{isLonLat}} \tab Test if an object has a longitude/latitude coordinate reference system \cr
    \code{\link{filename}} \tab Filename(s) to which a SpatRaster is linked \cr	
    \code{\link{compareGeom}} \tab Compare the geometry of SpatRasters \cr		
    \code{NAvalue} \tab Get or set the \code{NA} value (for reading from a file) \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}


\section{XI. Computing row, column, cell numbers and coordinates}{

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
    \code{\link{cellFromRowColCombine}} \tab Cell numbers from all combinations of row and column numbers \cr    \code{cellsFromExtent} \tab Cell numbers from extent object \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}

\section{XII. Writing files}{

  \tabular{ll}{
   \bold{Basic}\cr
    \code{\link{values}} \tab Assign new values to the cells of a SpatRaster \cr
    \code{\link{writeRaster}} \tab Write all values of SpatRaster to disk \cr
   .\cr
   \bold{Advanced}\cr
   \code{\link{blockSize}} \tab Get suggested block size for reading and writing \cr
   \code{\link{writeStart}} \tab Open a file for writing \cr
   \code{\link{writeValues}} \tab Write some values  \cr
   \code{\link{writeStop}} \tab Close the file after writing \cr
    --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}	


\section{XIII. Manipulation of SpatVector objects}{
The name in \bold{bold} is the equivalent command in ArcGIS. 

\tabular{ll}{
  \code{c} \tab \bold{append} combine "rbind" vector data of the same (vector) type \cr
  \code{erase} or "-" \tab \bold{erase} parts of a polygons\cr
  \code{intersect} or "*"  \tab \bold{intersect} polygons\cr
  \code{union} or "+" \tab \bold{union} polygons\cr
  \code{cover} \tab \bold{update} and \bold{identity} a polygons\cr
  \code{symdif} \tab\bold{symmetrical difference} of two polygons \cr 
  \code{aggregate} \tab \bold{dissolve} smaller polygons into larger ones \cr
  \code{disaggregate} \tab \bold{explode}: turn polygon parts into separate polygons \cr
  \code{crop} \tab \bold{clip} vector data using a rectangle (SpatExtent )\cr
  \code{select} \tab \bold{select} - interactively select spatial features\cr
  \code{click} \tab \bold{identify} attributes by clicking on a map\cr
  \code{merge} \tab \bold{Join table} \cr
  \code{extract} \tab spatial queries between SpatVector and SpatRaster objects \cr
  \code{as.data.frame} \tab coerce coordinates into a data.frame\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr  
}
}


\section{XIV. SpatExtent objects}{

  \tabular{ll}{
    \code{\link{extent}} \tab Create a SpatExtent object \cr
    \code{\link{intersect}} \tab Intersect two SpatExtent objects  \cr
    \code{\link{union}} \tab Combine two SpatExtent objects \cr
    \code{\link{Math-methods}} \tab round/floor/ceiling of the coordinates of a SpatExtent \cr			
    \code{\link{align}} \tab Align a SpatExtent with a SpatRaster \cr
    \code{\link{draw}} \tab Create a SpatExtent object by drawing it on top of a map (plot) \cr		
--------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}

\section{XV. Miscellaneous}{

\tabular{ll}{
    \code{\link{terraOptions}} \tab Show, set, save or get session options \cr
    \code{\link{tmpFiles}} \tab Show or remove temporary files \cr
    \code{\link{canProcessInMemory}} \tab Test whether a file can be created in memory \cr
    \code{\link{readStart}} \tab Open file connections for efficient multi-chunck reading \cr
    \code{\link{readStop}} \tab Close file connections \cr
    \code{\link{inMemory}} \tab Are the cell values in memory? \cr
    \code{\link{fromDisk}} \tab Are the cell values read from a file? \cr	
--------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}

\author{Except where indicated otherwise, the methods and functions in this package were written by Robert J. Hijmans.}

\section{Acknowledgments}{
This package stands on the shoulders of giants (notably for GDAL, GEOS, GeographicLib)
}

\keyword{ package }
\keyword{ spatial }


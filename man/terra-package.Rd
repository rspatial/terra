\name{terra-package}
\alias{terra-package}

\alias{terra}

\docType{package}

\title{terra}

\description{
This package is conceived as a replacement of the raster package. It has a very similar, but simpler, interface, and it is faster.

Like the raster package, "terra" provides classes and functions to manipulate geographic (spatial) data in "raster" format. Raster data divides space into cells (rectangles; pixels) of equal size (in units of the coordinate reference system). Such continuous spatial data are also referred to as "grid" data, and be contrasted with discrete (object based) spatial data such as points, lines, polygons (for which there is also support). An important difference is that "terra" has a single main class, "SpatRaster", as opposed to the three classes in the "raster" package (RasterLayer, RasterStack, RasterBrick). 

The package should be particularly useful when using very large datasets that can not be loaded into the computer's memory. Functions will work correctly, because they they process large files in chunks, i.e., they read, compute, and write blocks of data, without loading all values into memory at once. To speed up computations, all native computations are done in C++, exposed via a Rcpp module. Methods that heavily interact with native R functions (e.g. "predict") are written in R, calling C++ functions for reading/writing.


Below is a list of some of the most important functions grouped by theme.
}

\details{
 The package implements classes for spatial data (see \link{SpatRaster-class}) and supports
  \itemize{
    \item Creation of SpatRaster objects from scratch or from file
    \item Handling extremely large raster files
    \item Raster algebra and overlay functions
    \item Distance, neighborhood (focal) and patch functions
    \item Polygon, line and point to raster conversion
    \item Model predictions
    \item Summarizing raster values
    \item Easy access to raster cell-values
    \item Plotting (making maps)
    \item Manipulation of raster extent, resolution and origin
    \item Computation of row, col and cell numbers to coordinates and vice versa
    \item Reading and writing various raster file types
  }
. 
}


\section{I. Creating SpatRaster objects}{

SpatRasterobjects can be created, from scratch, files, or from objects of other classes, with the following functions:
  \tabular{ll}{
    \code{\link{rast}}\tab Create a SpatRaster from scratch, file, or another object\cr
    \code{\link{c}} \tab Combine SpatRasters (multiple layers) (like raster::stack)\cr
    \code{\link{subset}} \tab Select layers of a SpatRaster\cr
  --------------------------- \tab --------------------------------------------------------------------------------------------------- \cr  
}
}


\section{II. Changing the spatial extent and/or resolution of a SpatRaster}{

  \tabular{ll}{
    \code{\link{merge}} \tab Combine SpatRasters with different extents (but same origin and resolution) \cr
    \code{\link{mosaic}} \tab Combine SpatRasters with different extents and a function for the values in overlapping areas \cr
    \code{\link{crop}} \tab Select a geographic subset of a SpatRaster \cr
    \code{\link{extend}} \tab Enlarge a SpatRaster \cr
    \code{\link{trim}} \tab Trim a SpatRaster by removing exterior rows and/or columns that only have NAs\cr
    \code{\link{aggregate}} \tab Combine cells of a SpatRaster to create larger cells \cr
    \code{\link{disaggregate}} \tab Subdivide cells \cr
    \code{\link{resample}} \tab Warp values to a SpatRaster with a different origin or resolution \cr	
    \code{\link{warp}} \tab project values to a SpatRaster with a different coordinate reference system \cr	
    \code{\link{shift}} \tab Adjust the location of SpatRaster \cr
    \code{\link{flip}} \tab  Flip values horizontally or vertically \cr
    \code{\link{rotate}} \tab Rotate values around the date-line (for lon/lat data) \cr
    \code{\link{t}} \tab Transpose a SpatRaster\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}


\section{III. Raster algebra}{
	\tabular{ll}{
    \code{\link{Arith-methods}} \tab Arith functions (\code{+, -, *, ^, \%\%, \%/\%, /}) \cr
	\code{\link{Math-methods}} \tab Math functions like \code{abs, sqrt, trunc, log, log10, exp, sin, round} \cr	
    \code{\link{Logic-methods}} \tab Logic functions (\code{!, &, |}) \cr
    \code{\link{Summary-methods}} \tab Summary functions (\code{mean, max, min, range, prod, sum, any, all}) \cr
    \code{\link{Compare-methods}} \tab Compare functions (\code{==, !=, >, <, <=, >=}) \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
	}
}  
  
  
\section{IV. Cell based computation}{

  \tabular{ll}{
    \code{\link{reduce}} \tab Computations on a single SpatRaster (similar to "calc" in the "raster" package) \cr
    \code{\link{overlay}} \tab Computations on multiple SpatRaster objects \cr
    \code{\link{cover}} \tab First layer covers second layer except where the first layer is \code{NA} \cr
    \code{\link{mask}} \tab Use values from first SpatRaster except where cells of the mask SpatRaster are \code{NA} (or another value of choice)\cr
    \code{\link{cut}} \tab Reclassify values using ranges \cr
    \code{\link{subs}} \tab Reclassify values using an 'is-becomes' matrix \cr	
    \code{\link{reclassify}} \tab Reclassify using a 'from-to-becomes' matrix \cr
    \code{\link{init}} \tab Initialize cells with new values \cr	
    \code{\link{rastApply}} \tab Computations on groups of layers in SpatRaster \cr
    \code{\link{rastSelect}} \tab Select cell values from different layers using an index SpatRaster\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}


\section{V. Spatial contextual computation}{

  \tabular{ll}{
    \code{\link{distance}} \tab Shortest distance to a cell that is not \code{NA}\cr
    \code{\link{gridDistance}} \tab Distance when traversing grid cells that are not \code{NA} \cr
    \code{\link{distanceFromPoints}} \tab Shortest distance to any point in a set of points \cr
    \code{\link{direction}} \tab Direction (azimuth) to or from cells that are not \code{NA}\cr
    \code{\link{focal}} \tab Focal (neighborhood; moving window) functions \cr
    \code{\link{localFun}} \tab Local association (using neighborhoods) functions \cr
    \code{\link{boundaries}} \tab Detection of boundaries (edges)\cr
    \code{\link{clump}} \tab Find clumps (patches) \cr
    \code{\link{adjacent}} \tab Identify cells that are adjacent to a set of cells of a SpatRaster \cr
    \code{\link{area}} \tab Compute area of cells (for longitude/latitude data) \cr
    \code{\link{terrain}} \tab Compute slope, aspect and other characteristics from elevation data \cr	
    \code{\link{Moran}} \tab Compute global or local Moran or Geary indices of spatial autocorrelation  \cr	
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

You can coerce SpatRasters to Spatial* objects using \code{as}, as in \code{as(object, 'SpatialGridDataFrame')}  

  \tabular{ll}{
    \code{\link{rast}} \tab SpatRaster from matrix and other objects\cr
    \code{\link{rasterize}} \tab Rasterizing points, lines or polygons\cr
    \code{\link{rasterToPoints}} \tab Create points from a SpatRaster \cr
    \code{\link{rasterToPolygons}} \tab Create polygons from a SpatRaster \cr
    \code{\link{rasterToContour}} \tab Contour lines from a SpatRaster \cr
    \code{\link{rasterFromXYZ}} \tab Createa a SpatRaster from regularly spaces points\cr
    \code{\link{rasterFromCells}} \tab Create a SpatRaster from anoter Raster object and cell numbers\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}


\section{VIII. Summarizing}{

  \tabular{ll}{
    \code{\link{cellStats}} \tab Summarize SpatRaster cell values with a function \cr
    \code{\link{summary}} \tab Summary of the values of a SpatRaster (quartiles and mean) \cr
    \code{\link{freq}} \tab Frequency table of SpatRaster cell values \cr
    \code{\link{crosstab}} \tab Cross-tabulate two SpatRasters\cr
    \code{\link{unique}} \tab Get the unique values in a SpatRaster \cr
    \code{\link{zonal}} \tab Summarize a SpatRaster by zones in another SpatRaster \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}


\section{IX. Accessing values of SpatRaster cells}{

Apart from the function listed below, you can also use indexing with \code{[} for cell numbers, and \code{[[} for row / column number combinations \cr

\tabular{ll}{
    \code{\link{values}} \tab Get all cell values (fails with very large rasters), or a row of values (safer) \cr
    \code{\link{as.matrix}} \tab Get cell values as a matrix \cr
    \code{\link{as.array}} \tab Get cell values as an array \cr
    \code{\link{extract}} \tab Extract cell values from a SpatRaster (e.g., by cell, coordinates, polygon)\cr
    \code{\link{sampleRandom}} \tab Random sample \cr
    \code{\link{sampleRegular}} \tab Regular sample \cr
    \code{\link{minValue}} \tab Get the minimum value of the cells of a SpatRaster (not always known) \cr
    \code{\link{maxValue}} \tab Get the maximum value of the cells of a SpatRaster (not always known) \cr
    \code{\link{setMinMax}} \tab Compute the minimum and maximum value of a SpatRaster if these are not known \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}


\section{X. Plotting}{

  \tabular{ll}{
    \bold{Maps} \tab \cr
    \code{\link{plot}} \tab Plot a SpatRaster. The main method to create a map \cr
    \code{\link{plotRGB}} \tab Combine three layers (red, green, blue channels) into a single 'real color' image \cr
    \code{\link{image}} \tab Plot a SpatRaster with the image function \cr
    \code{\link{persp}} \tab Perspective plot of a SpatRaster \cr
    \code{\link{contour}} \tab Contour plot of a SpatRaster \cr
    \code{\link{filledContour}} \tab Filled contour plot of a SpatRaster \cr
    \code{\link{text}} \tab Plot the values of a SpatRaster on top of a map \cr
    .\cr
    \bold{Interacting with a map} \tab \cr
    \code{\link{zoom}} \tab Zoom in to a part of a map \cr
    \code{\link{click}} \tab Query values of SpatRaster or Spatial* objects by clicking on a map \cr
    \code{\link{select}} \tab Select a geometric subset of a SpatRaster or Spatial* object \cr
    \code{\link{drawPoly}} \tab Create a SpatialPolygons object by drawing it \cr
    \code{\link{drawLine}} \tab Create a SpatialLines object by drawing it \cr
    \code{\link{drawExtent}} \tab Create an Extent object by drawing it \cr	
    .\cr
	\bold{Other plots} \tab \cr
    \code{\link{plot}} \tab x-y scatter plot of the values of two SpatRaster objects\cr
    \code{\link{hist}} \tab Histogram of SpatRaster values \cr
    \code{\link{barplot}} \tab barplot of a SpatRaster \cr
    \code{\link{density}} \tab Density plot of SpatRaster values \cr
    \code{\link{pairs}} \tab Pairs plot for layers in a SpatRaster \cr
    \code{\link{boxplot}} \tab Box plot of the values of one or multiple layers\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}

\section{XI. Getting and setting SpatRaster dimensions }{

Basic parameters of existing SpatRasters can be obtained, and in most cases changed. If there are values associated with a SpatRaster object (either in memory or via a link to a file) these are lost when you change the number of columns or rows or the resolution. This is not the case when the extent is changed (as the number of columns and rows will not be affected). Similarly, with \bold{projection} you can set the projection, but this does not transform the data (see \link{warp} for that).   

  \tabular{ll}{
    \code{\link{ncol}}\tab The number of columns \cr
    \code{\link{nrow}} \tab The number of rows \cr
    \code{\link{ncell}} \tab The number of cells (can not be set directly, only via ncol or nrow) \cr
    \code{\link{res}} \tab The resolution (x and y) \cr
    \code{\link{nlyr}} \tab How many layers does the object have? \cr	
    \code{\link{names}} \tab Get or set the layer names \cr	
    \code{\link{xres}} \tab The x resolution (can be set with res) \cr
    \code{\link{yres}} \tab The y resolution (can be set with res)\cr
    \code{\link{xmin}} \tab The minimum x coordinate (or longitude) \cr
    \code{\link{xmax}} \tab The maximum x coordinate (or longitude) \cr
    \code{\link{ymin}} \tab The minimum y coordinate (or latitude) \cr
    \code{\link{ymax}} \tab The maximum y coordinate (or latitude) \cr
	\code{\link{extent}} \tab The extent (minimum and maximum x and y coordinates) \cr
	\code{\link{origin}} \tab The origin of a SpatRaster\cr
    \code{\link{crs}} \tab The coordinate reference system (map projection) \cr	
    \code{\link{isLonLat}} \tab Test if an object has a longitude/latitude coordinate reference system \cr
    \code{\link{filename}} \tab Filename(s) to which a SpatRaster is linked \cr	
    \code{\link{compareRaster}} \tab Compare the geometry of SpatRasters \cr		
    \code{\link{NAvalue}} \tab Get or set the \code{NA} value (for reading from a file) \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
	
}
}


\section{XII. Computing row, column, cell numbers and coordinates}{

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
    \code{\link{cellsFromExtent}} \tab Cell numbers from extent object \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}

\section{XIII. Writing files}{

  \tabular{ll}{
   \bold{Basic}\cr
    \code{\link{setValues}} \tab Put new values in a SpatRaster \cr
    \code{\link{writeRaster}} \tab Write all values of SpatRaster to disk \cr
    \code{\link{KML}} \tab Save raster as KML file \cr	
   .\cr
   \bold{Advanced}\cr
   \code{\link{blockSize}} \tab Get suggested block size for reading and writing \cr
   \code{\link{writeStart}} \tab Open a file for writing \cr
   \code{\link{writeValues}} \tab Write some values  \cr
   \code{\link{writeStop}} \tab Close the file after writing \cr
    --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}	

\section{XV. SpatExtent objects}{

  \tabular{ll}{
    \code{\link{extent}} \tab Create an extent object \cr
    \code{\link{intersect}} \tab Intersect two extent objects  \cr
    \code{\link{union}} \tab Combine two extent objects \cr
    \code{\link{round}} \tab round/floor/ceiling of the coordinates of an Extent object \cr			
    \code{\link{alignExtent}} \tab Align an extent with a SpatRaster \cr
    \code{\link{drawExtent}} \tab Create an Extent object by drawing it on top of a map (see plot) \cr		
    --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}

\section{XVI. Miscellaneous}{

\tabular{ll}{
    \code{\link{terraOptions}} \tab Show, set, save or get session options \cr
    \code{\link{pointDistance}} \tab Distance between points \cr
    \code{\link{hdr}} \tab Write header file for a number of raster formats \cr
    \code{\link{extension}} \tab Get or set the extension of a filename \cr
    \code{\link{cv}} \tab Coefficient of variation  \cr
    \code{\link{modal}} \tab Modal value  \cr
    \code{\link{sampleInt}} \tab Random sample of (possibly very large) range of integer values \cr
    \code{\link{showTmpFiles}} \tab Show temporary files \cr
    \code{\link{removeTmpFiles}} \tab Remove temporary files \cr
    --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}

\section{XVII. For programmers}{

\tabular{ll}{
  \code{\link{canProcessInMemory}} \tab Test whether a file can be created in memory \cr
  \code{\link{readStart}} \tab Open file connections for efficient multi-chunck reading \cr
  \code{\link{readStop}} \tab Close file connections \cr
  \code{\link{rasterTmpFile}}  \tab Get a name for a temporary file \cr
  \code{\link{inMemory}} \tab Are the cell values in memory? \cr
  \code{\link{fromDisk}} \tab Are the cell values read from a file? \cr	
    --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}

\author{Except where indicated otherwise, the methods and functions in this package were written by Robert J. Hijmans}

\section{Acknowledgments}{
to do
}



\keyword{ package }
\keyword{ spatial }


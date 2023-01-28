\name{terra-package}
\alias{terra-package}

\alias{terra}

\docType{package}

\title{Description of the methods in the terra package}

\description{

\code{terra} provides methods to manipulate geographic (spatial) data in "raster" and "vector" form. Raster data divide space into rectangular grid cells and they are commonly used to represent spatially continuous phenomena, such as elevation or the weather. Satellite images also have this data structure, and in that context grid cells are often referred to as pixels. In contrast, "vector" spatial data (points, lines, polygons) are typically used to represent discrete spatial entities, such as a road, country, or bus stop. 

The package implements two main classes (data types): \code{SpatRaster} and \code{SpatVector}. \code{SpatRaster} supports handling large raster files that cannot be loaded into memory; local, focal, zonal, and global raster operations; polygon, line and point to raster conversion; integration with modeling methods to make spatial predictions; and more. \code{SpatVector} supports all types of geometric operations such as intersections.

Additional classes include \code{SpatExtent}, which is used to define a spatial extent (bounding box); \code{SpatRasterDataset}, which represents a collection of sub-datasets for the same area. Each sub-dataset is a SpatRaster with possibly many layers, and may, for example, represent different weather variables; and \code{SpatRasterCollection} and \code{SpatVectorCollection} that are equivalent to lists of \code{SpatRaster} or \code{SpatVector} objects. There is also a \code{SpatGraticule} class to assist in adding a longitude/latitude lines and labels to a map with another coordinate reference system.

These classes hold a C++ pointer to the data "reference class" and that creates some limitations. They cannot be recovered from a saved R session either or directly passed to nodes on a computer cluster. Generally, you should use \code{\link{writeRaster}} to save \code{SpatRaster} objects to disk (and pass a filename or cell values to cluster nodes). Also see \code{\link{wrap}}. 

The \code{terra} package is conceived as a replacement of the \code{raster} package. \code{terra} has a very similar, but simpler, interface, and it is faster than \code{raster}. At the bottom of this page there is a table that shows differences in the methods between the two packages. 

Below is a list of some of the most important methods grouped by theme. 

---------------------------------------------------------------------------------------------------------------------
}


\section{\bold{SpatRaster}}{

---------------------------------------------------------------------------------------------------------------------}

\section{I. Creating, combining and sub-setting}{

  \tabular{ll}{
    \code{\link{rast}}\tab Create a SpatRaster from scratch, file, or another object\cr
    \code{\link{c}} \tab Combine SpatRasters (multiple layers)\cr
    \code{\link{add<-}} \tab Add a SpatRaster to another one\cr
    \code{\link{subset}} or \code{[[}, or \code{$} \tab Select layers of a SpatRaster\cr
    \code{\link{selectRange}} \tab Select cell values from different layers using an index layer\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}


\section{II. Changing the spatial extent or resolution}{
Also see the methods in section VIII

  \tabular{ll}{
    \code{\link{merge}} \tab Combine SpatRasters with different extents (but same origin and resolution) \cr
    \code{\link{mosaic}} \tab Combine SpatRasters with different extents using a function for overlapping cells \cr
    \code{\link{crop}} \tab Select a geographic subset of a SpatRaster \cr
    \code{\link{extend}} \tab Add rows and/or columns to a SpatRaster \cr
    \code{\link{trim}} \tab Trim a SpatRaster by removing exterior rows and/or columns that only have NAs\cr
    \code{\link{aggregate}} \tab Combine cells of a SpatRaster to create larger cells \cr
    \code{\link{disagg}} \tab Subdivide cells \cr
    \code{\link{resample}} \tab Resample (warp) values to a SpatRaster with a different origin and/or resolution \cr 
	\code{\link{project}} \tab Project (warp) values to a SpatRaster with a different coordinate reference system \cr
    \code{\link{shift}} \tab Adjust the location of SpatRaster \cr
    \code{\link{flip}} \tab  Flip values horizontally or vertically \cr
    \code{\link{rotate}} \tab Rotate values around the date-line (for lon/lat data) \cr
    \code{\link{t}} \tab Transpose a SpatRaster\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}

  
\section{III. Local (cell based) methods}{

\subsection{Apply-like methods}{

  \tabular{ll}{
    \code{\link{app}} \tab Apply a function to all cells, across layers, typically to summarize (as in \code{base::apply}) \cr
    \code{\link{tapp}} \tab Apply a function to groups of layers (as in \code{base::tapply} and \code{stats::aggregate})\cr
    \code{\link{lapp}} \tab Apply a function to using the layers of a SpatRaster as variables\cr
    \code{\link{sapp}} \tab Apply a function to each layer\cr
    \code{\link{rapp}} \tab Apply a function to a spatially variable range of layers\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
	}
}

\subsection{Arithmetic, logical, and standard math methods}{

  \tabular{ll}{
    \code{\link{Arith-methods}} \tab Standard arithmetic methods (\code{+, -, *, ^, \%\%, \%/\%, /}) \cr
    \code{\link{Compare-methods}} \tab Comparison methods for SpatRaster (\code{==, !=, >, <, <=, >=m is.na, is.finite}) \cr 	
    \code{\link{not.na}} \tab a one-step equivalent to \code{!is.na} \cr
    \code{\link{Summary-methods}} \tab \code{mean, max, min, median, sum, range, prod,} \cr
    \tab \code{any, all, stdev, which.min, which.max, anyNA, noNA, allNA}\cr
    \code{\link{Logic-methods}} \tab Boolean methods (\code{!, &, |}) \cr
	\code{\link{Math-methods}} \tab \code{abs, sign, sqrt, ceiling, floor, trunc, cummax, cummin, cumprod,} \cr
	\tab \code{cumsum, log, log10, log2, log1p, acos, acosh, asin, asinh, atan, atanh,} \cr	 
	\tab \code{exp, expm1, cos, cosh, sin, sinh, tan, tanh, round, signif}\cr
    \code{\link{as.bool}}\tab create a Boolean (logical) SpatRaster \cr
    \code{\link{as.int}}\tab create an integer (whole numbers) SpatRaster \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
	}
}

\subsection{Other methods}{

  \tabular{ll}{
    \code{\link{approximate}} \tab Compute missing values for cells by interpolation across layers \cr
    \code{\link{roll}} \tab Rolling functions such as the rolling mean \cr
    \code{\link{cellSize}} \tab Compute the area of cells \cr
    \code{\link{classify}} \tab (Re-)classify values \cr
    \code{\link{subst}} \tab Substitute (replace) cell values \cr
    \code{\link{cover}} \tab First layer covers second layer except where the first layer is \code{NA} \cr
    \code{\link{init}} \tab Initialize cells with new values \cr
    \code{\link{mask}} \tab Replace values in a SpatRaster based on values in another SpatRaster\cr
    \code{\link{which.lyr}} \tab which is the first layer that is \code{TRUE}?\cr
    \code{\link{segregate}} \tab Make a 0/1 layer for each unique value \cr
    \code{\link{rangeFill}} \tab Make a 0/1 SpatRaster for a time series \cr
	
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
  }
}

\section{IV. Zonal and global methods}{
  \tabular{ll}{
    \code{\link{expanse}} \tab Compute the summed area of cells \cr
    \code{\link{crosstab}} \tab Cross-tabulate two SpatRasters\cr
    \code{\link{freq}} \tab Frequency table of SpatRaster cell values \cr
    \code{\link{global}} \tab Summarize SpatRaster cell values with a function \cr
    \code{\link{quantile}} \tab Quantiles \cr
    \code{\link{layerCor}} \tab Correlation between layers \cr
    \code{\link{stretch}} \tab Stretch values \cr
    \code{\link{scale}} \tab Scale values \cr
    \code{\link[terra]{summary}} \tab Summary of the values of a SpatRaster (quartiles and mean) \cr
    \code{\link{unique}} \tab Get the unique values in a SpatRaster \cr
    \code{\link{zonal}} \tab Summarize a SpatRaster by zones in another SpatRaster \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}


\section{V. Situation (spatial context) based methods}{

  \tabular{ll}{
    \code{\link{adjacent}} \tab Identify cells that are adjacent to a set of cells of a SpatRaster \cr
    \code{\link{boundaries}} \tab Detection of boundaries (edges)\cr
    \code{\link{distance}} \tab Shortest distance to a cell that is not \code{NA} or to or from a vector object\cr
    \code{\link{gridDist}} \tab Shortest distance through adjacent grid cells\cr
    \code{\link{costDist}} \tab Shortest distance considering cell-varying friction \cr
    \code{\link{direction}} \tab Direction (azimuth) to or from cells that are not \code{NA}\cr
    \code{\link{focal}} \tab Focal (neighborhood; moving window) functions \cr
    \code{\link{focal3D}} \tab Three dimensional (row, col, lyr) focal functions \cr
    \code{\link{focalCpp}} \tab Faster focal by using custom C++ functions \cr
    \code{\link{focalReg}} \tab Regression between layers for focal areas \cr
    \code{\link{focalPairs}} \tab Apply a function (e.g. a correlation coeffient) to focal values for pairs of layers \cr
    \code{\link{patches}} \tab Find patches (clumps) \cr
    \code{\link{sieve}} \tab Sieve filter to remove small patches\cr
    \code{\link{terrain}} \tab Compute slope, aspect and other terrain characteristics from elevation data \cr
    \code{\link{viewshed}} \tab Compute vieshed (showing areas that are visible from a particular location \cr
    \code{\link{shade}} \tab Compute hill shade from slope and aspect layers \cr
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
    \code{\link{values}} \tab cell values (fails with very large rasters)\cr
    \code{\link{values<-}} \tab Set new values to the cells of a SpatRaster \cr
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


\section{VIII. Getting and setting dimensions }{

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
	\code{\link{ext}} \tab Get or set the extent (minimum and maximum x and y coordinates ("bounding box") \cr
	\code{\link{origin}} \tab The origin of a SpatRaster\cr
    \code{\link{crs}} \tab The coordinate reference system (map projection) \cr
    \code{\link{is.lonlat}} \tab Test if an object has (or may have) a longitude/latitude coordinate reference system\cr
    \code{\link{sources}} \tab Get the filename(s) to which a SpatRaster is linked \cr
    \code{\link{inMemory}} \tab Are the data sources in memory (or on disk)? \cr
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


\section{X. Time related methods}{

\tabular{ll}{
    \code{\link{time}} \tab Get or set time\cr
    \code{\link{fillTime}} \tab can add empty layers in between existing layers to assure that the time step between layers is constant \cr
    \code{\link{mergeTime}} \tab combine multiple rasters, perhaps partly overlapping in time, into a single time series\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}


\section{XI. Methods for categorical rasters}{

\tabular{ll}{
    \code{\link{is.factor}} \tab Are there categorical layers?\cr
    \code{\link{levels}} \tab Get active categories, or set categories\cr
    \code{\link{activeCat}} \tab Get or set the active category\cr
    \code{\link{cats}} \tab Get categories (active and inactive)\cr
    \code{\link{set.cats}} \tab Set categories in place \cr
    \code{\link{concats}} \tab Combine SpatRasters with different categories\cr
    \code{\link{catalyze}} \tab Create a layer for each category \cr
    \code{\link{as.numeric}} \tab use the active category to create a non-categorical SpatRaster\cr
    \code{\link{as.factor}} \tab Make the layers of a SpatRaster categorical\cr


  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}


\section{XII. Writing SpatRaster files}{

\subsection{Basic}{
  \tabular{ll}{
    \code{\link{writeRaster}} \tab Write all values of SpatRaster to disk. You can set the filetype, datatype, compression. \cr
    \code{\link{writeCDF}} \tab Write SpatRaster data to a netCDF file\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
 	}
}

\subsection{Advanced}{
  \tabular{ll}{
   \code{\link{readStart}} \tab Open file connections for efficient multi-chunk reading \cr
   \code{\link{readStop}} \tab Close file connections \cr
   \code{\link{writeStart}} \tab Open a file for writing \cr
   \code{\link{writeValues}} \tab Write some values  \cr
   \code{\link{writeStop}} \tab Close the file after writing \cr
   \code{\link{blocks}} \tab Get blocksize for reading files (when not writing) \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
\cr
}
}
}

\section{XIII. Miscellaneous SpatRaster methods}{

\tabular{ll}{
    \code{\link{terraOptions}} \tab Show, set, or get session options, mostly to control memory use and to set write options\cr
    \code{\link{sources}} \tab Show the data sources of a SpatRaster \cr
    \code{\link{tmpFiles}} \tab Show or remove temporary files \cr
    \code{\link{mem_info}} \tab memory needs and availability \cr
    \code{\link{inMemory}} \tab Are the cell values in memory? \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}


\section{XIV. SpatRasterDataset}{

A SpatRasterDataset contains SpatRaster objects that are sub-datasets for the same area. They all have the same extent and resolution.

\tabular{ll}{
  \code{\link{sds}} \tab Create a SpatRasterDataset from a file with subdatasets (ncdf or hdf) or from SpatRaster objects \cr
  \code{[} or \code{$} \tab Extract a SpatRaster \cr
  \code{\link{names}} \tab Get the names of the sub-datasets \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}

\section{XV. SpatRasterCollections}{
A SpatRasterCollection is a vector of SpatRaster objects. Unlike for a SpatRasterDataset, there the extent and resolution of the SpatRasters do not need to match each other. 

\tabular{ll}{
\code{\link{sprc}} \tab create a SpatRasterCollection from (a list of) SpatRaster objects\cr
\code{\link{length}} \tab how many SpatRasters does the SpatRasterCollection have?\cr
\code{\link{crop}} \tab crop a SpatRasterCollection\cr
\code{\link{impose}} \tab force the members of SpatRasterCollection to the same geometry\cr
\code{\link{merge}} \tab merge the members of a SpatRasterCollection\cr
\code{\link{mosaic}} \tab mosaic (merge with a function for overlapping areas) the members of a SpatRasterCollection\cr
\code{\link{[}} \tab extract a SpatRaster\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}}



\section{\bold{SpatVector}}{

---------------------------------------------------------------------------------------------------------------------}

\section{XVI. Create SpatVector objects}{
\tabular{ll}{
  \code{\link{vect}} \tab Create a SpatVector from a file (for example a "shapefile") or from another object\cr
  \code{\link{vector_layers}} \tab list or delete layers in a vector database such as GPGK\cr
  \code{rbind} \tab append SpatVectors of the same geometry type\cr
  \code{\link{unique}} \tab remove duplicates \cr
  \code{\link{na.omit}} \tab remove empty geometries and/or fields that are \code{NA} \cr
  \code{\link{project}} \tab Project a SpatVector to a different coordinate reference system \cr
  \code{\link{writeVector}} \tab Write SpatVector data to disk \cr
  \code{\link{centroids}} \tab Get the centroids of a SpatVector\cr
  \code{\link{voronoi}} \tab Voronoi diagram \cr
  \code{\link{delaunay}} \tab Delaunay triangles\cr
  \code{\link{convHull}} \tab Compute the convex hull of a SpatVector \cr  
  \code{\link{fillHoles}}\tab Remove or extract holes from polygons\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr  
}}


\section{XVII. Properties of SpatVector objects}{
\tabular{ll}{
  \code{\link{geom}} \tab returns the geometries as matrix or WKT\cr
  \code{\link{crds}} \tab returns the coordinates as a matrix\cr
  \code{\link{linearUnits}} \tab returns the linear units of the crs (in meter)\cr
  \code{\link{ncol}}\tab The number of columns (of the attributes)\cr
  \code{\link{nrow}} \tab The number of rows (of the geometries and attributes)\cr
  \code{\link{names}} \tab Get or set the layer names \cr
  \code{\link{ext}} \tab Get the extent (minimum and maximum x and y coordinates ("bounding box") \cr
  \code{\link{crs}} \tab The coordinate reference system (map projection) \cr
  \code{\link{is.lonlat}} \tab Test if an object has (or may have) a longitude/latitude coordinate reference system\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr  
}}


\section{XVIII. Geometric queries}{
\tabular{ll}{
  \code{\link{adjacent}} \tab find adjacent polygons\cr
  \code{\link{expanse}} \tab computes the area covered by polygons\cr
  \code{\link{nearby}} \tab find nearby geometries\cr
  \code{\link{nearest}} \tab find the nearest geometries\cr
  \code{\link{relate}} \tab geometric relationships such as "intersects", "overlaps", and "touches"\cr
  \code{\link{perim}} \tab computes the length of the perimeter of polygons, and the length of lines\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr  
}}


\section{XIX. Geometric operations}{

\tabular{ll}{
  \code{\link{erase}} or "-" \tab erase (parts of) geometries\cr
  \code{\link{intersect}} or "*"  \tab intersect geometries\cr
  \code{\link{union}} or "+" \tab Merge geometries\cr
  \code{\link{cover}} \tab update polygons\cr
  \code{\link{symdif}} \tab symmetrical difference of two polygons \cr 
  \code{\link{aggregate}} \tab dissolve smaller polygons into larger ones \cr
  \code{\link{buffer}} \tab buffer geometries \cr
  \code{\link{disagg}} \tab split multi-geometries into separate geometries \cr
  \code{\link{crop}} \tab clip geometries using a rectangle (SpatExtent) or SpatVector \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr  
}}


\section{XX. SpatVector attributes}{
We use the term "attributes" for the tabular data (data.frame) associated with vector geometries.

\tabular{ll}{
  \code{\link{extract}} \tab spatial queries between SpatVector and SpatVector (e.g. point in polygons) \cr
  \code{\link{sel}} \tab select - interactively select geometries\cr
  \code{\link{click}} \tab identify attributes by clicking on a map\cr
  \code{\link{merge}} \tab Join a table with a SpatVector \cr
  \code{\link{as.data.frame}} \tab get attributes as a data.frame\cr
  \code{\link{as.list}} \tab get attributes as a list\cr
  \code{\link{values}} \tab Get the attributes of a SpatVector \cr
  \code{\link{values<-}} \tab Set new attributes to the geometries of a SpatRaster \cr  
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr  
}}


\section{XXI. Change geometries (for display, experimentation)}{

\tabular{ll}{
  \code{\link{shift}} \tab change the position geometries by shifting their coordinates in horizontal and/or vertical direction\cr
  \code{\link{spin}}\tab rotate geometries around an origin\cr
  \code{\link{rescale}} \tab shrink (or expand) geometries, for example to make an inset map \cr
  \code{\link{flip}} \tab flip geometries vertically or horizontally\cr
  \code{\link{t}} \tab transpose geometries (switch x and y)\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr  
}}


\section{XXII. Geometry properties and topology}{

\tabular{ll}{
  \code{\link{width}} \tab the minimum diameter of the geometries \cr
  \code{\link{clearance}} \tab the minimum clearance of the geometries \cr
  \code{\link{sharedPaths}} \tab shared paths (arcs) between line or polygon geometries\cr
  \code{\link{simplifyGeom}} \tab simplify geometries\cr
  \code{\link{gaps}} \tab find gaps between polygon geometries \cr
  \code{\link{fillHoles}} \tab get or remove the polygon holes\cr
  \code{\link{makeNodes}} \tab create nodes on lines \cr
  \code{\link{mergeLines}} \tab connect lines to form polygons  \cr
  \code{\link{removeDupNodes}} \tab remove duplicate nodes in geometries and optionally rounds the coordinates  \cr
  \code{\link{is.valid}} \tab check if geometries are valid \cr
  \code{\link{makeValid}} \tab attempt to repair invalid geometries \cr  
  \code{\link{snap}} \tab make boundaries of geometries identical if they are very close to each other  \cr
  \code{\link{erase} (single argument)} \tab remove parts of geometries that overlap \cr
  \code{\link{union} (single argument)} \tab create new polygons such that there are no overlapping polygons \cr
  \code{\link{rotate}} \tab rotate to (dis-) connect them across the date-line \cr
  \code{\link{normalize.longitude}} \tab move geomtries that are outside of the -180 to 180 degrees range.  \cr
  \code{\link{elongate}} \tab make lines longer by extending both sides \cr
  \code{\link{combineGeoms}} \tab combine geometries that overlap, share a border, or are within a minimum distance of each other \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}}


\section{XXIII. SpatVectorCollections}{
A SpatVectorCollection is a vector of SpatVector objects. 

\tabular{ll}{
\code{\link{svc}} \tab create a SpatVectorCollection from (a list of) SpatVector objects\cr
\code{\link{length}} \tab how many SpatRasters does the SpatRasterCollection have?\cr
\code{\link{[}} \tab extract a SpatVector\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}}


\section{\bold{Other classes}}{

---------------------------------------------------------------------------------------------------------------------}


\section{XXIV. SpatExtent}{
  \tabular{ll}{
    \code{\link{ext}} \tab Create a SpatExtent object. For example to \code{\link{crop}} a Spatial dataset\cr
    \code{\link{intersect}} \tab Intersect two SpatExtent objects, same as \code{-}  \cr
    \code{\link{union}} \tab Combine two SpatExtent objects, same as \code{+} \cr
    \code{\link{Math-methods}} \tab round/floor/ceiling of a SpatExtent \cr
    \code{\link{align}} \tab Align a SpatExtent with a SpatRaster \cr
    \code{\link{draw}} \tab Create a SpatExtent by drawing it on top of a map (plot) \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}

\section{XXV. SpatGraticule}{
  \tabular{ll}{
    \code{\link{graticule}} \tab Create a graticule\cr
    \code{\link{crop}} \tab crop a graticule\cr
	\code{\link[=plot,SpatGraticule,missing-method]{plot<SpatGraticule>}} \tab plot a graticule \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}


\section{\bold{General methods}}{

---------------------------------------------------------------------------------------------------------------------}

\section{XXVI. Conversion between spatial data objects from different packages}{

You can coerce SpatRasters to Raster* objects, after loading the \code{raster} package, with \code{as(object, "Raster")}, or \code{raster(object)} or \code{brick(object)} or \code{stack(object)}
 
  \tabular{ll}{
    \code{\link{rast}} \tab SpatRaster from matrix and other objects\cr
    \code{\link{vect}} \tab SpatVector from \code{sf} or \code{Spatial*} vector data\cr
    \code{sf::st_as_sf} \tab sf object from SpatVector\cr
    \code{\link{rasterize}} \tab Rasterizing points, lines or polygons\cr
    \code{\link{rasterizeWin}} \tab Rasterize points with a moving window\cr
    \code{\link{rasterizeGeom}} \tab Rasterize attributes of geometries such as "count", "area", or "length"\cr
    \code{\link{as.points}} \tab Create points from a SpatRaster or SpatVector \cr
    \code{\link{as.lines}} \tab Create points from a SpatRaster or SpatVector\cr
    \code{\link{as.polygons}} \tab Create polygons from a SpatRaster \cr
    \code{\link{as.contour}}  \tab Contour lines from a SpatRaster \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
  }
}

\section{XXVII. Plotting}{

\subsection{Maps}{
  \tabular{ll}{
    \code{\link{plot}} \tab Plot a SpatRaster or SpatVector. The main method to create a map \cr
	\code{\link{points}} \tab Add points to a map\cr
	\code{\link{lines}} \tab Add lines to a map\cr
	\code{\link{polys}} \tab Add polygons to a map \cr
    \code{\link{text}} \tab Add text (such as the values of a SpatRaster or SpatVector) to a map \cr
    \code{\link{image}} \tab Alternative to plot to make a map with a SpatRaster \cr
    \code{\link{plotRGB}} \tab Combine three layers (red, green, blue channels) into a single "real color" plot \cr
	\code{\link[=plot,SpatGraticule,missing-method]{plot<SpatGraticule>}} \tab plot a graticule \cr
    \code{\link{sbar}} \tab Add a scalebar to a map\cr
    \code{\link{north}} \tab Add a north arrow to a map\cr
    \code{\link{inset}} \tab Add a small inset (overview) map\cr
    \code{\link{dots}} \tab Make a dot-density map \cr
    \code{\link{cartogram}} \tab Make a cartogram \cr
    \code{\link{persp}} \tab Perspective plot of a SpatRaster \cr
    \code{\link{contour}} \tab Contour plot or filled-contour plot of a SpatRaster \cr
    \code{\link{colorize}} \tab Combine three layers (red, green, blue channels) into a single layer with a color-table \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
	}
}

\subsection{Interacting with a map}{
  \tabular{ll}{
    \code{\link{zoom}} \tab Zoom in to a part of a map by drawing a bounding box on it \cr
    \code{\link{click}} \tab Query values of SpatRaster or SpatVector by clicking on a map \cr
    \code{\link{sel}} \tab Select a spatial subset of a SpatRaster or SpatVector by drawing on a map\cr
    \code{\link{draw}} \tab Create a SpatExtent or SpatVector by drawing on a map \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
	}
}

\subsection{Other plots}{
  \tabular{ll}{
    \code{\link{plot}} \tab x-y scatter plot of the values of (a sample of) the layers of two SpatRaster objects\cr
    \code{\link{hist}} \tab Histogram of SpatRaster values \cr
    \code{\link{barplot}} \tab Bar plot of a SpatRaster \cr
    \code{\link{density}} \tab Density plot of SpatRaster values \cr
    \code{\link{pairs}} \tab Pairs plot for layers in a SpatRaster \cr
    \code{\link{boxplot}} \tab Box plot of the values of a SpatRaster\cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
	}
}
}

\section{\bold{Comparison with the raster package}}{

---------------------------------------------------------------------------------------------------------------------}

\section{XXVIII. New method names}{

\code{terra} has a single class \code{SpatRaster} for which \code{raster} has three (\code{RasterLayer, RasterStack, RasterBrick}). Likewise there is a single class for vector data \code{SpatVector} that replaces six \code{Spatial*} classes. Most method names are the same, but note the following important differences in methods names with the \code{raster} package

\tabular{ll}{
\bold{raster package} \tab \bold{terra package}\cr
\code{raster, brick, stack}\tab\code{\link{rast}}\cr
\code{rasterFromXYZ}\tab\code{\link{rast}( , type="xyz")}\cr
\code{stack, addLayer}\tab\code{\link{c}}\cr
\code{addLayer} \tab\code{\link{add<-}}\cr
\code{area} \tab\code{\link{cellSize}} or \code{\link{expanse}}\cr
\code{approxNA} \tab\code{\link{approximate}}\cr
\code{calc}\tab\code{\link{app}}\cr
\code{cellFromLine, cellFromPolygon,}\tab \code{\link{cells}}\cr
\code{cellsFromExtent}\tab\code{\link{cells}}\cr
\code{cellStats}\tab\code{\link{global}}\cr
\code{clump}\tab\code{\link{patches}}\cr
\code{compareRaster}\tab\code{\link{compareGeom}}\cr
\code{corLocal}\tab\code{\link{focalCor}}\cr
\code{coordinates}\tab\code{\link{crds}}\cr
\code{couldBeLonLat}\tab\code{\link{is.lonlat}}\cr
\code{disaggregate} \tab\code{\link{disagg}}\cr
\code{distanceFromPoints} \tab \code{\link{distance}}\cr
\code{drawExtent, drawPoly, drawLine} \tab \code{\link{draw}}\cr
\code{dropLayer}\tab\code{\link{subset}}\cr
\code{extent}\tab\code{\link{ext}}\cr
\code{getValues}\tab\code{\link{values}}\cr
\code{isLonLat, isGlobalLonLat}\tab\code{\link{is.lonlat}}\cr
\code{layerize}\tab\code{\link{segregate}}\cr
\code{layerStats}\tab\code{\link{layerCor}}\cr
\code{movingFun}\tab\code{\link{roll}}\cr
\code{NAvalue}\tab\code{\link{NAflag}}\cr
\code{nlayers}\tab\code{\link{nlyr}}\cr
\code{overlay}\tab\code{\link{lapp}}\cr
\code{projectRaster}\tab\code{\link{project}}\cr
\code{rasterToPoints}\tab\code{\link{as.points}}\cr
\code{rasterToPolygons}\tab\code{\link{as.polygons}}\cr
\code{reclassify, subs, cut}\tab\code{\link{classify}}\cr
\code{sampleRandom, sampleRegular}\tab\code{\link{spatSample}}\cr
\code{shapefile}\tab\code{\link{vect}}\cr
\code{stackApply}\tab\code{\link{tapp}}\cr
\code{stackSelect}\tab\code{\link{selectRange}}\cr
}
}

\section{XXIX. Changed behavior}{

Also note that even if function names are the same in \code{terra} and \code{raster}, their output can be different. In most cases this was done to get more consistency in the returned values (and thus fewer errors in the downstream code that uses them). It other cases it simply seemed better. Here are some examples:

\tabular{ll}{
\code{\link{as.polygons}}\tab By default, \code{terra} returns dissolved polygons\cr
\code{\link{quantile}}\tab computes by cell, across layers instead of the other way around\cr
\code{\link{extract}}\tab By default, \code{terra} returns a matrix, with the first column the sequential ID of the vectors. \cr
\tab \code{raster} returns a list (for lines or polygons) or a matrix (for points, but without the ID\cr
\tab column. You can use \code{list=TRUE} to get the results as a list\cr
\code{\link{values}}\tab \code{terra} always returns a matrix. \code{raster} returns a vector for a \code{RasterLayer}\cr
\code{\link{Summary-methods}}\tab With \code{raster}, \code{mean(x, y)} and \code{mean(stack(x, y)} return the same result, a single\cr
\tab layer with the mean of all cell values. This is also what \code{terra} returns with \cr
\tab \code{mean(c(x, y))}, but with \code{mean(x, y)} the parallel mean is returned -- that is, the\cr
\tab computation is done layer-wise, and the number of layers in the output is the same as\cr
\tab  that of \code{x} and \code{y} (or the larger of the two if they are not the same). This affects \cr
\tab all summary functions (\code{sum}, \code{mean}, \code{median}, \code{which.min}, \code{which.max}, \code{min}, \code{max},\cr
\tab \code{prod}, \code{any}, \code{all}, \code{stdev}), except \code{range}, which is not implemented for this case \cr
\tab (you can use \code{min} and \code{max} instead) \cr
  --------------------------- \tab ------------------------------------------------------------------------------------------ \cr
}
}

\section{Authors}{
Except where indicated otherwise, the methods and functions in this package were written by Robert Hijmans. The configuration scripts were written by Roger Bivand. Some of code using the GEOS library was adapted from code by Edzer Pebesma for \code{sf}. Michael Sumner contributed various bits and pieces.
}

\section{Acknowledgments}{
This package is an attempt to climb on the shoulders of giants (GDAL, PROJ, GEOS, NCDF, GeographicLib, Rcpp, R). Many people have contributed by asking questions or \href{https://github.com/rspatial/terra}{raising issues}. Feedback and suggestions by Márcia Barbosa, Kendon Bell, Andrew Gene Brown, Jean-Luc Dupouey, Krzysztof Dyba, Alex Ilich, Gerald Nelson, Jakub Nowosad and Monika Tomaszewska have been especially helpful.
}


\keyword{ package }
\keyword{ spatial }


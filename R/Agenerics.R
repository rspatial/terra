#if (!isGeneric("#")) { setGeneric("#", function(object) standardGeneric("#")) }


if (!isGeneric("cols")) {setGeneric("cols", function(x, ...) standardGeneric("cols"))}
if (!isGeneric("window<-")) {setGeneric("window<-", function(x, ..., value) standardGeneric("window<-"))}

# terra
if (!isGeneric("meta")) { setGeneric("meta", function(x, ...) standardGeneric("meta"))}	

if (!isGeneric("app")) { setGeneric("app", function(x, ...) standardGeneric("app"))}
if (!isGeneric("lapp")) { setGeneric("lapp", function(x, ...) standardGeneric("lapp"))}
if (!isGeneric("rapp")) { setGeneric("rapp", function(x, ...) standardGeneric("rapp"))}
if (!isGeneric("tapp")) { setGeneric("tapp", function(x, ...) standardGeneric("tapp"))}

if (!isGeneric("add<-")) {setGeneric("add<-", function(x, value) standardGeneric("add<-"))}
if (!isGeneric("align")) { setGeneric("align", function(x, y, ...) standardGeneric("align"))}
if (!isGeneric("arith")) {setGeneric("arith", function(x,...) standardGeneric("arith"))}
if (!isGeneric("as.contour")) {setGeneric("as.contour", function(x,...) standardGeneric("as.contour"))}	
if (!isGeneric("as.lines")) {setGeneric("as.lines", function(x,...) standardGeneric("as.lines"))}	
if (!isGeneric("as.points")) {setGeneric("as.points", function(x,...) standardGeneric("as.points"))}	
if (!isGeneric("as.polygons")) {setGeneric("as.polygons", function(x,...) standardGeneric("as.polygons"))}	
if (!isGeneric("classify")) { setGeneric("classify", function(x, ...) standardGeneric("classify")) }
if (!isGeneric("cells")) { setGeneric("cells", function(x, y, ...) standardGeneric("cells")) }
if (!isGeneric("collapse")) {setGeneric("collapse", function(x, ...) standardGeneric("collapse"))}
if (!isGeneric("compareGeom")) {setGeneric("compareGeom", function(x,y,...) standardGeneric("compareGeom"))}	
if (!isGeneric("crosstab")) { setGeneric("crosstab", function(x, y, ...) standardGeneric("crosstab")) }

if (!isGeneric("depth")) {setGeneric("depth", function(x,...) standardGeneric("depth"))}
if (!isGeneric("depth<-")) {setGeneric("depth<-", function(x, value) standardGeneric("depth<-"))}
if (!isGeneric("draw")) {setGeneric("draw", function(x,...) standardGeneric("draw"))}
if (!isGeneric("ext")) { setGeneric("ext", function(x, ...)	standardGeneric("ext"))}
if (!isGeneric("ext<-")) { setGeneric("ext<-", function(x, value) standardGeneric("ext<-")) }	
if (!isGeneric("fill") ) { setGeneric("fill", function(x, ...) standardGeneric("fill")) }
if (!isGeneric("geomtype")) {setGeneric("geomtype", function(x,...) standardGeneric("geomtype"))}

if (!isGeneric("global")) {setGeneric("global", function(x, ...) standardGeneric("global"))}	
if (!isGeneric("is.points")) {setGeneric("is.points", function(x,...) standardGeneric("is.points"))}
if (!isGeneric("is.lines")) {setGeneric("is.lines", function(x,...) standardGeneric("is.lines"))}
if (!isGeneric("is.polygons")) {setGeneric("is.polygons", function(x,...) standardGeneric("is.polygons"))}

if (!isGeneric("isTRUE")) { setGeneric("isTRUE", function(x) standardGeneric("isTRUE"))}
if (!isGeneric("isFALSE")) { setGeneric("isFALSE", function(x) standardGeneric("isFALSE"))}
if (!isGeneric("minmax")) {setGeneric("minmax", function(x) standardGeneric("minmax"))}	
if (!isGeneric("nsrc")) { setGeneric("nsrc", function(x) standardGeneric("nsrc")) }	
if (!isGeneric("perimeter")) {setGeneric("perimeter", function(x, ...) standardGeneric("perimeter"))}	
if (!isGeneric("project")) {setGeneric("project", function(x,...) standardGeneric("project"))}	
if (!isGeneric("pack")) {setGeneric("pack", function(x, ...) standardGeneric("pack"))}	
if (!isGeneric("rats")) { setGeneric("rats", function(x, ...) standardGeneric("rats")) }
if (!isGeneric("rast") ) { setGeneric("rast", function(x, ...) standardGeneric("rast")) }
if (!isGeneric("rev") ) { setGeneric("rev", function(x) standardGeneric("rev")) }
if (!isGeneric("sds") ) { setGeneric("sds", function(x, ...) standardGeneric("sds")) }
if (!isGeneric("select")) {setGeneric("select", function(x, ...) standardGeneric("select"))}	
if (!isGeneric("separate")) {setGeneric("separate", function(x, ...) standardGeneric("separate"))}	

if (!isGeneric("selectRange")) {setGeneric("selectRange", function(x, ...) standardGeneric("selectRange"))}
if (!isGeneric("setValues")) {setGeneric("setValues", function(x, ...) standardGeneric("setValues"))}	
if (!isGeneric("slope")) {setGeneric("slope", function(x, ...) standardGeneric("slope"))}	
if (!isGeneric("size")) {setGeneric("size", function(x, ...) standardGeneric("size"))}	
if (!isGeneric("sources")) {setGeneric("sources", function(x, ...) standardGeneric("sources"))}	
if (!isGeneric("spatSample")) { setGeneric("spatSample", function(x, ...) standardGeneric("spatSample"))}
if (!isGeneric("time")) {setGeneric("time", function(x,...) standardGeneric("time"))}
if (!isGeneric("time<-")) {setGeneric("time<-", function(x, value) standardGeneric("time<-"))}
if (!isGeneric("nlyr")) { setGeneric("nlyr", function(x) standardGeneric("nlyr")) }	
if (!isGeneric("nlyr<-")) { setGeneric("nlyr<-", function(x, ..., value) standardGeneric("nlyr<-")) }	
if (!isGeneric("vect") ) { setGeneric("vect", function(x, ...) standardGeneric("vect")) }


## shared with "raster"
if (!isGeneric("%in%")) { setGeneric("%in%", function(x, table)	standardGeneric("%in%")) }	
if (!isGeneric("adjacent")) {setGeneric("adjacent", function(x, ...) standardGeneric("adjacent"))}	
if (!isGeneric("animate")) { setGeneric("animate", function(x, ...) standardGeneric("animate")) }	

if (!isGeneric("area")) {setGeneric("area", function(x, ...) standardGeneric("area"))}	
if (!isGeneric("as.data.frame")) { setGeneric("as.data.frame", function(x, ...) standardGeneric("as.data.frame"))}	

if (!isGeneric("as.list")) { setGeneric("as.list", function(x, ...) standardGeneric("as.list"))}	
if (!isGeneric("as.factor")) {setGeneric("as.factor", function(x) standardGeneric("as.factor"))}
if (!isGeneric("atan2")) { setGeneric("atan2", function(y, x) standardGeneric("atan2"))}	
if (!isGeneric("barplot")) {setGeneric("barplot", function(height,...) standardGeneric("barplot"))}	
if (!isGeneric("bbox")) {setGeneric("bbox", function(obj) standardGeneric("bbox"))}	
if (!isGeneric("boundaries")) {	setGeneric("boundaries", function(x, ...) standardGeneric("boundaries"))}	
if (!isGeneric("boxplot")) { setGeneric("boxplot", function(x, ...) standardGeneric("boxplot"))}
if (!isGeneric("buffer")) {setGeneric("buffer", function(x, ...) standardGeneric("buffer"))}	
if (!isGeneric("clamp")) { setGeneric("clamp", function(x, ...) standardGeneric("clamp")) }
if (!isGeneric("click")) {setGeneric("click", function(x, ...)standardGeneric("click"))}	
if (!isGeneric("contour")) { setGeneric("contour", function(x,...) standardGeneric("contour"))}	
if (!isGeneric("cover")) {setGeneric("cover", function(x, y, ...) standardGeneric("cover"))}			
if (!isGeneric("crop")) { setGeneric("crop", function(x, y, ...) standardGeneric("crop")) }	
if (!isGeneric("crs")) { setGeneric("crs", function(x, ...)	standardGeneric("crs")) }	
if (!isGeneric("crs<-")) { setGeneric("crs<-", function(x, ..., value) standardGeneric("crs<-")) }	
if (!isGeneric("density")) { setGeneric("density", function(x, ...) standardGeneric("density"))}	
if (!isGeneric("distance")) {setGeneric("distance", function(x, y, ...)standardGeneric("distance"))}	
if (!isGeneric("extract")) { setGeneric("extract", function(x, y, ...) standardGeneric("extract"))}	
if (!isGeneric("expand")) {setGeneric("expand", function(x, y, ...) standardGeneric("expand"))}	
if (!isGeneric("flip")) {setGeneric("flip", function(x, ...) standardGeneric("flip")) }
if (!isGeneric("focal")) { setGeneric("focal", function(x, ...) standardGeneric("focal")) }
if (!isGeneric("freq")) { setGeneric("freq", function(x, ...) standardGeneric("freq")) }	
if (!isGeneric("geom")) { setGeneric("geom", function(x,...) standardGeneric("geom"))}	
#if (!isGeneric("gridDistance")) {setGeneric("gridDistance", function(x, ...) standardGeneric("gridDistance"))}	
if (!isGeneric("hasValues")) {setGeneric("hasValues", function(x, ...) standardGeneric("hasValues")) }	
if (!isGeneric("head")) { setGeneric("head", function(x, ...) standardGeneric("head"))}	
if (!isGeneric("ifel")) {setGeneric("ifel", function(test, yes, no, ...) standardGeneric("ifel"))}	
if (!isGeneric("image")) {setGeneric("image", function(x, ...)standardGeneric("image"))}	
if (!isGeneric("init")) {setGeneric("init", function(x, ...) standardGeneric("init"))}	
if (!isGeneric("interpolate")) { setGeneric("interpolate", function(object, ...) standardGeneric("interpolate"))}	
if (!isGeneric("is.factor")) {setGeneric("is.factor", function(x) standardGeneric("is.factor"))}
if (!isGeneric("isLonLat")) { setGeneric("isLonLat", function(x, ...) standardGeneric("isLonLat"))}	
if (!isGeneric("mask")) { setGeneric("mask", function(x, mask, ...) standardGeneric("mask")) }	
if (!isGeneric("match")) { setGeneric("match", function(x, table, nomatch=NA_integer_, incomparables=NULL)		standardGeneric("match"))}
if (!isGeneric("modal")) {setGeneric("modal", function(x, ...) standardGeneric("modal"))}
if (!isGeneric("mosaic")) {setGeneric("mosaic", function(x, ...) standardGeneric("mosaic"))}

if (!isGeneric("ncell")) { setGeneric("ncell", function(x) standardGeneric("ncell")) }	
if (!isGeneric("nrow")) { setGeneric("nrow", function(x) standardGeneric("nrow")) }	
if (!isGeneric("ncol")) { setGeneric("nrow", function(x) standardGeneric("nrow")) }	
if (!isGeneric("ncol<-")) { setGeneric("ncol<-", function(x, ..., value) standardGeneric("ncol<-")) }	
if (!isGeneric("nrow<-")) { setGeneric("nrow<-", function(x, ..., value) standardGeneric("nrow<-")) }	
if (!isGeneric("origin")) {	setGeneric("origin", function(x, ...) standardGeneric("origin")) }
if (!isGeneric("pairs")) { setGeneric("pairs", function(x, ...)	standardGeneric("pairs"))}
if (!isGeneric("persp")) { setGeneric("persp", function(x,...) standardGeneric("persp")) }	
if (!isGeneric("plot")) { setGeneric("plot", function(x, y,...) standardGeneric("plot"))}	
if (!isGeneric("plotRGB")) { setGeneric("plotRGB", function(x, ...)standardGeneric("plotRGB"))}	
if (!isGeneric("predict")) {setGeneric("predict", function(object, ...) standardGeneric("predict"))}
if (!isGeneric("quantile")) {setGeneric("quantile", function(x, ...)standardGeneric("quantile"))}	
if (!isGeneric("rasterize")) {setGeneric("rasterize", function(x, y, ...) standardGeneric("rasterize"))}
if (!isGeneric("readStart")) {setGeneric("readStart", function(x, ...) standardGeneric("readStart"))}
if (!isGeneric("readStop")) {setGeneric("readStop", function(x)	standardGeneric("readStop"))}
if (!isGeneric("res")) { setGeneric("res", function(x) standardGeneric("res")) }
if (!isGeneric("res<-")) { setGeneric("res<-", function(x, value) standardGeneric("res<-")) }
if (!isGeneric("rectify")) {setGeneric("rectify", function(x, ...) standardGeneric("rectify"))}	
if (!isGeneric("resample")) { setGeneric("resample", function(x, y, ...) standardGeneric("resample"))}
if (!isGeneric("select")) {setGeneric("select", function(x, ...) standardGeneric("select"))}
if (!isGeneric("setMinMax")) {setGeneric("setMinMax", function(x) standardGeneric("setMinMax"))}	
if (!isGeneric("shift")) {setGeneric("shift", function(x, ...) standardGeneric("shift"))}	
if (!isGeneric("stdev")) { setGeneric("stdev", function(x, ...) standardGeneric("stdev")) }	
if (!isGeneric("subset")) {setGeneric("subset", function(x, ...) standardGeneric("subset")) }
if (!isGeneric("summary")) {setGeneric("summary", function(object, ...) standardGeneric("summary")) }
if (!isGeneric("t")) { setGeneric("t", function(x) standardGeneric("t"))}	
if (!isGeneric("tail")) { setGeneric("tail", function(x, ...) standardGeneric("tail"))}	
if (!isGeneric("text")) { setGeneric("text", function(x, ...) standardGeneric("text")) }	
if (!isGeneric("transpose")) { setGeneric("transpose", function(x, ...) standardGeneric("transpose"))}	
if (!isGeneric("trim")) { setGeneric("trim", function(x, ...) standardGeneric("trim")) }	
if (!isGeneric("xres")) { setGeneric("xres", function(x) standardGeneric("xres")) }
if (!isGeneric("yres")) { setGeneric("yres", function(x) standardGeneric("yres")) }
if (!isGeneric("zonal")) {setGeneric("zonal", function(x, z, ...) standardGeneric("zonal"))}	

if (!isGeneric("yFromRow")) { setGeneric("yFromRow", function(object, row) standardGeneric("yFromRow")) }
if (!isGeneric("xFromCol")) { setGeneric("xFromCol", function(object, col) standardGeneric("xFromCol")) }                  
if (!isGeneric("colFromX")) { setGeneric("colFromX", function(object, x) standardGeneric("colFromX")) }                  
if (!isGeneric("rowFromY")) { setGeneric("rowFromY", function(object, y) standardGeneric("rowFromY")) }                  
if (!isGeneric("cellFromXY")) { setGeneric("cellFromXY", function(object, xy) standardGeneric("cellFromXY")) }
if (!isGeneric("cellFromRowCol")) { setGeneric("cellFromRowCol", function(object, row, col, ...) standardGeneric("cellFromRowCol")) }
if (!isGeneric("cellFromRowColCombine")) { setGeneric("cellFromRowColCombine", function(object, row, col, ...) standardGeneric("cellFromRowColCombine")) }

if (!isGeneric("xyFromCell")) { setGeneric("xyFromCell", function(object, cell, ...) standardGeneric("xyFromCell")) } 
if (!isGeneric("yFromCell")) { setGeneric("yFromCell", function(object, cell) standardGeneric("yFromCell")) }             
if (!isGeneric("xFromCell")) { setGeneric("xFromCell", function(object, cell) standardGeneric("xFromCell")) }               
if (!isGeneric("rowColFromCell")) { setGeneric("rowColFromCell", function(object, cell) standardGeneric("rowColFromCell")) }
if (!isGeneric("rowFromCell")) { setGeneric("rowFromCell", function(object, cell) standardGeneric("rowFromCell")) } 
if (!isGeneric("colFromCell")) { setGeneric("colFromCell", function(object, cell) standardGeneric("colFromCell")) }

if (!isGeneric("readStart")) { setGeneric("readStart", function(x, ...) standardGeneric("readStart")) }	
if (!isGeneric("readValues")) { setGeneric("readValues", function(x, ...) standardGeneric("readValues")) }	
if (!isGeneric("readStop")) { setGeneric("readStop", function(x, ...) standardGeneric("readStop")) }


if (!isGeneric("setMinMax")) {setGeneric("setMinMax", function(x, ...) standardGeneric("setMinMax")) }	
if (!isGeneric("stretch")) {setGeneric("stretch", function(x, ...) standardGeneric("stretch")) }	

if (!isGeneric("unique")) { setGeneric("unique", function(x, incomparables=FALSE, ...) standardGeneric("unique")) }	
if (!isGeneric("values")) { setGeneric("values", function(x, ...) standardGeneric("values")) }	
if (!isGeneric("values<-")) { setGeneric("values<-", function(x, value) standardGeneric("values<-"))}	
if (!isGeneric("which.max")) {setGeneric("which.max", function(x) standardGeneric("which.max"))}	
if (!isGeneric("which.min")) {setGeneric("which.min", function(x) standardGeneric("which.min"))}	
#if (!isGeneric("Which")) {setGeneric("Which", function(x) standardGeneric("Which"))}	

if (!isGeneric("writeStart")) {	setGeneric("writeStart", function(x, filename, ...)	standardGeneric("writeStart")) }
if (!isGeneric("writeStop")) { setGeneric("writeStop", function(x, ...) standardGeneric("writeStop")) }
#if (!isGeneric("writeVals")) { setGeneric("writeVals", function(x, v, ...) standardGeneric("writeVals")) }
if (!isGeneric("writeValues")) { setGeneric("writeValues", function(x, v, ...) standardGeneric("writeValues")) }


if (!isGeneric("writeRaster")) {setGeneric("writeRaster", function(x, filename, ...) standardGeneric("writeRaster"))}
if (!isGeneric("writeVector")) {setGeneric("writeVector", function(x, filename, ...) standardGeneric("writeVector"))}
if (!isGeneric("xmin")) {setGeneric("xmin", function(x) standardGeneric("xmin"))}
if (!isGeneric("xmax")) {setGeneric("xmax", function(x)	standardGeneric("xmax"))}
if (!isGeneric("ymin")) {setGeneric("ymin", function(x)	standardGeneric("ymin"))}
if (!isGeneric("ymax")) {setGeneric("ymax", function(x)	standardGeneric("ymax"))}
if (!isGeneric("xmin<-")) { setGeneric("xmin<-", function(x, ..., value) standardGeneric("xmin<-"))}	
if (!isGeneric("xmax<-")) { setGeneric("xmax<-", function(x, ..., value) standardGeneric("xmax<-"))}	
if (!isGeneric("ymin<-")) { setGeneric("ymin<-", function(x, ..., value) standardGeneric("ymin<-"))}	
if (!isGeneric("ymax<-")) { setGeneric("ymax<-", function(x, ..., value) standardGeneric("ymax<-"))}	




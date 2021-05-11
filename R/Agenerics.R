#if (!isGeneric("#")) { setGeneric("#", function(object) standardGeneric("#")) }
if (!isGeneric("weighted.mean")) {setGeneric("weighted.mean", function(x, w, ...) standardGeneric("weighted.mean"))}

if (!isGeneric("split")) {setGeneric("split", function(x, ...) standardGeneric("split"))}
if (!isGeneric("cellSize")) {setGeneric("cellSize", function(x, ...) standardGeneric("cellSize"))}

if (!isGeneric("sharedPaths")) {setGeneric("sharedPaths", function(x, ...) standardGeneric("sharedPaths"))}
if (!isGeneric("isTRUE")) {setGeneric("isTRUE", function(x) standardGeneric("isTRUE"))}
if (!isGeneric("subst")) {setGeneric("subst", function(x, ...) standardGeneric("subst"))}
if (!isGeneric("RGB2col")) {setGeneric("RGB2col", function(x, ...) standardGeneric("RGB2col"))}
if (!isGeneric("RGB")) {setGeneric("RGB", function(x, ...) standardGeneric("RGB"))}
if (!isGeneric("RGB<-")) {setGeneric("RGB<-", function(x, ..., value) standardGeneric("RGB<-"))}
if (!isGeneric("autocor")) {setGeneric("autocor", function(x, ...) standardGeneric("autocor"))}
if (!isGeneric("delauny")) {setGeneric("delauny", function(x, ...) standardGeneric("delauny"))}
if (!isGeneric("voronoi")) {setGeneric("voronoi", function(x, ...) standardGeneric("voronoi"))}
if (!isGeneric("convHull")) {setGeneric("convHull", function(x, ...) standardGeneric("convHull"))}
if (!isGeneric("relate")) {setGeneric("relate", function(x, y, ...) standardGeneric("relate"))}
if (!isGeneric("intersect")) {setGeneric("intersect", function(x, y, ...) standardGeneric("intersect"))}
if (!isGeneric("erase")) {setGeneric("erase", function(x, y, ...) standardGeneric("erase"))}
if (!isGeneric("nearby")) {setGeneric("nearby", function(x, ...) standardGeneric("nearby"))}
if (!isGeneric("nearest")) {setGeneric("nearest", function(x, ...) standardGeneric("nearest"))}
if (!isGeneric("cartogram")) {setGeneric("cartogram", function(x, ...) standardGeneric("cartogram"))}
if (!isGeneric("dots")) {setGeneric("dots", function(x, ...) standardGeneric("dots"))}
if (!isGeneric("crds")) {setGeneric("crds", function(x, ...) standardGeneric("crds"))}
if (!isGeneric("symdif")) {setGeneric("symdif", function(x, y, ...) standardGeneric("symdif"))}
if (!isGeneric("union")) {setGeneric("union", function(x, y, ...) standardGeneric("union"))}
if (!isGeneric("median")) {setGeneric("median", function(x, na.rm) standardGeneric("median"))}
if (!isGeneric("polys")) {setGeneric("polys", function(x,...) standardGeneric("polys"))}
if (!isGeneric("centroids")) {setGeneric("centroids", function(x, ...) standardGeneric("centroids"))}
if (!isGeneric("coltab")) {setGeneric("coltab", function(x, ...) standardGeneric("coltab"))}
if (!isGeneric("coltab<-")) { setGeneric("coltab<-", function(x, ..., value) standardGeneric("coltab<-")) }
if (!isGeneric("copy")) { setGeneric("copy", function(x, ...) standardGeneric("copy")) }
if (!isGeneric("window")) {setGeneric("window", function(x, ...) standardGeneric("window"))}
if (!isGeneric("window<-")) {setGeneric("window<-", function(x, ..., value) standardGeneric("window<-"))}
if (!isGeneric("NAflag")) {setGeneric("NAflag", function(x, ...) standardGeneric("NAflag"))}
if (!isGeneric("NAflag<-")) {setGeneric("NAflag<-", function(x, ..., value) standardGeneric("NAflag<-"))}
if (!isGeneric("app")) { setGeneric("app", function(x, ...) standardGeneric("app"))}
if (!isGeneric("lapp")) { setGeneric("lapp", function(x, ...) standardGeneric("lapp"))}
if (!isGeneric("rapp")) { setGeneric("rapp", function(x, ...) standardGeneric("rapp"))}
if (!isGeneric("tapp")) { setGeneric("tapp", function(x, ...) standardGeneric("tapp"))}
if (!isGeneric("sapp")) { setGeneric("sapp", function(x, ...) standardGeneric("sapp"))}
if (!isGeneric("add<-")) {setGeneric("add<-", function(x, value) standardGeneric("add<-"))}
if (!isGeneric("align")) { setGeneric("align", function(x, y, ...) standardGeneric("align"))}
if (!isGeneric("as.contour")) {setGeneric("as.contour", function(x,...) standardGeneric("as.contour"))}
if (!isGeneric("as.lines")) {setGeneric("as.lines", function(x,...) standardGeneric("as.lines"))}
if (!isGeneric("as.points")) {setGeneric("as.points", function(x,...) standardGeneric("as.points"))}
if (!isGeneric("as.polygons")) {setGeneric("as.polygons", function(x,...) standardGeneric("as.polygons"))}
if (!isGeneric("classify")) { setGeneric("classify", function(x, ...) standardGeneric("classify")) }
if (!isGeneric("cells")) { setGeneric("cells", function(x, y, ...) standardGeneric("cells")) }
if (!isGeneric("tighten")) {setGeneric("tighten", function(x, ...) standardGeneric("tighten"))}
if (!isGeneric("compareGeom")) {setGeneric("compareGeom", function(x,y,...) standardGeneric("compareGeom"))}
if (!isGeneric("crosstab")) { setGeneric("crosstab", function(x, y, ...) standardGeneric("crosstab")) }
if (!isGeneric("describe")) { setGeneric("describe", function(x, ...) standardGeneric("describe"))}
if (!isGeneric("depth")) {setGeneric("depth", function(x,...) standardGeneric("depth"))}
if (!isGeneric("depth<-")) {setGeneric("depth<-", function(x, value) standardGeneric("depth<-"))}
if (!isGeneric("draw")) {setGeneric("draw", function(x,...) standardGeneric("draw"))}
if (!isGeneric("ext")) { setGeneric("ext", function(x, ...)	standardGeneric("ext"))}
if (!isGeneric("ext<-")) { setGeneric("ext<-", function(x, value) standardGeneric("ext<-")) }
if (!isGeneric("fillHoles") ) { setGeneric("fillHoles", function(x, ...) standardGeneric("fillHoles")) }
if (!isGeneric("geomtype")) {setGeneric("geomtype", function(x, ...) standardGeneric("geomtype"))}
if (!isGeneric("datatype")) {setGeneric("datatype", function(x, ...) standardGeneric("datatype"))}
if (!isGeneric("global")) {setGeneric("global", function(x, ...) standardGeneric("global"))}
if (!isGeneric("is.valid")) {setGeneric("is.valid", function(x,...) standardGeneric("is.valid"))}
if (!isGeneric("is.points")) {setGeneric("is.points", function(x,...) standardGeneric("is.points"))}
if (!isGeneric("is.lines")) {setGeneric("is.lines", function(x,...) standardGeneric("is.lines"))}
if (!isGeneric("is.polygons")) {setGeneric("is.polygons", function(x,...) standardGeneric("is.polygons"))}
if (!isGeneric("makeTiles")) {setGeneric("makeTiles", function(x,...) standardGeneric("makeTiles"))}
if (!isGeneric("vrt")) {setGeneric("vrt", function(x,...) standardGeneric("vrt"))}
if (!isGeneric("isTRUE")) { setGeneric("isTRUE", function(x) standardGeneric("isTRUE"))}
if (!isGeneric("isFALSE")) { setGeneric("isFALSE", function(x) standardGeneric("isFALSE"))}
if (!isGeneric("varnames")) {setGeneric("varnames", function(x,...) standardGeneric("varnames"))}
if (!isGeneric("varnames<-")) {setGeneric("varnames<-", function(x, value) standardGeneric("varnames<-"))}
if (!isGeneric("log")) {setGeneric("log", function(x,...) standardGeneric("log"))}
if (!isGeneric("longnames")) {setGeneric("longnames", function(x,...) standardGeneric("longnames"))}
if (!isGeneric("longnames<-")) {setGeneric("longnames<-", function(x, value) standardGeneric("longnames<-"))}
if (!isGeneric("minmax")) {setGeneric("minmax", function(x) standardGeneric("minmax"))}
if (!isGeneric("nsrc")) { setGeneric("nsrc", function(x) standardGeneric("nsrc")) }
if (!isGeneric("perim")) {setGeneric("perim", function(x, ...) standardGeneric("perim"))}
if (!isGeneric("project")) {setGeneric("project", function(x,...) standardGeneric("project"))}
if (!isGeneric("wrap")) {setGeneric("wrap", function(x, ...) standardGeneric("wrap"))}
if (!isGeneric("cats")) { setGeneric("cats", function(x, ...) standardGeneric("cats")) }
if (!isGeneric("setCats")) { setGeneric("setCats", function(x, ...) standardGeneric("setCats")) }
if (!isGeneric("as.raster")) { setGeneric("as.raster", function(x, ...) standardGeneric("as.raster"))}
if (!isGeneric("rast") ) { setGeneric("rast", function(x, ...) standardGeneric("rast")) }
if (!isGeneric("rev") ) { setGeneric("rev", function(x) standardGeneric("rev")) }
if (!isGeneric("sds") ) { setGeneric("sds", function(x, ...) standardGeneric("sds")) }
if (!isGeneric("svc") ) { setGeneric("svc", function(x, ...) standardGeneric("svc")) }
if (!isGeneric("sel")) {setGeneric("sel", function(x, ...) standardGeneric("sel"))}
if (!isGeneric("segregate")) {setGeneric("segregate", function(x, ...) standardGeneric("segregate"))}
if (!isGeneric("selectRange")) {setGeneric("selectRange", function(x, ...) standardGeneric("selectRange"))}
if (!isGeneric("setValues")) {setGeneric("setValues", function(x, ...) standardGeneric("setValues"))}
if (!isGeneric("expanse")) {setGeneric("expanse", function(x, ...) standardGeneric("expanse"))}
if (!isGeneric("size")) {setGeneric("size", function(x, ...) standardGeneric("size"))}
if (!isGeneric("sources")) {setGeneric("sources", function(x, ...) standardGeneric("sources"))}
if (!isGeneric("spatSample")) { setGeneric("spatSample", function(x, ...) standardGeneric("spatSample"))}
if (!isGeneric("terrain")) {setGeneric("terrain", function(x, ...) standardGeneric("terrain"))}
if (!isGeneric("time")) {setGeneric("time", function(x,...) standardGeneric("time"))}
if (!isGeneric("time<-")) {setGeneric("time<-", function(x, value) standardGeneric("time<-"))}
if (!isGeneric("nlyr")) { setGeneric("nlyr", function(x) standardGeneric("nlyr")) }
if (!isGeneric("nlyr<-")) { setGeneric("nlyr<-", function(x, ..., value) standardGeneric("nlyr<-")) }
if (!isGeneric("linearUnits")) {setGeneric("linearUnits", function(x, ...) standardGeneric("linearUnits"))}
if (!isGeneric("units")) {setGeneric("units", function(x) standardGeneric("units"))}
if (!isGeneric("units<-")) {setGeneric("units<-", function(x,value) standardGeneric("units<-"))}
if (!isGeneric("vect") ) { setGeneric("vect", function(x, ...) standardGeneric("vect")) }
if (!isGeneric("writeCDF")) {setGeneric("writeCDF", function(x, filename, ...) standardGeneric("writeCDF"))}
if (!isGeneric("writeVector")) {setGeneric("writeVector", function(x, filename, ...) standardGeneric("writeVector"))}
## shared with "raster"
if (!isGeneric("%in%")) { setGeneric("%in%", function(x, table)	standardGeneric("%in%")) }
if (!isGeneric("adjacent")) {setGeneric("adjacent", function(x, ...) standardGeneric("adjacent"))}
if (!isGeneric("animate")) { setGeneric("animate", function(x, ...) standardGeneric("animate")) }
if (!isGeneric("area")) {setGeneric("area", function(x, ...) standardGeneric("area"))}
if (!isGeneric("as.data.frame")) { setGeneric("as.data.frame", function(x, ...) standardGeneric("as.data.frame"))}
if (!isGeneric("as.list")) { setGeneric("as.list", function(x, ...) standardGeneric("as.list"))}
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
if (!isGeneric("extend")) {setGeneric("extend", function(x, y, ...) standardGeneric("extend"))}
if (!isGeneric("flip")) {setGeneric("flip", function(x, ...) standardGeneric("flip")) }
if (!isGeneric("focal")) { setGeneric("focal", function(x, ...) standardGeneric("focal")) }
if (!isGeneric("focalValues")) { setGeneric("focalValues", function(x, ...) standardGeneric("focalValues")) }
if (!isGeneric("freq")) { setGeneric("freq", function(x, ...) standardGeneric("freq")) }
if (!isGeneric("geom")) { setGeneric("geom", function(x,...) standardGeneric("geom"))}
if (!isGeneric("hasValues")) {setGeneric("hasValues", function(x, ...) standardGeneric("hasValues")) }
if (!isGeneric("head")) { setGeneric("head", function(x, ...) standardGeneric("head"))}
if (!isGeneric("ifel")) {setGeneric("ifel", function(test, ...) standardGeneric("ifel"))}
if (!isGeneric("image")) {setGeneric("image", function(x, ...)standardGeneric("image"))}
if (!isGeneric("init")) {setGeneric("init", function(x, ...) standardGeneric("init"))}
if (!isGeneric("inset")) {setGeneric("inset", function(x, ...) standardGeneric("inset"))}
if (!isGeneric("interpolate")) { setGeneric("interpolate", function(object, ...) standardGeneric("interpolate"))}
if (!isGeneric("is.factor")) {setGeneric("is.factor", function(x) standardGeneric("is.factor"))}
if (!isGeneric("is.lonlat")) { setGeneric("is.lonlat", function(x, ...) standardGeneric("is.lonlat"))}
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
if (!isGeneric("patches")) {setGeneric("patches", function(x, ...) standardGeneric("patches"))}
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
if (!isGeneric("spin")) {setGeneric("spin", function(x, ...) standardGeneric("spin"))}
if (!isGeneric("rotate")) {setGeneric("rotate", function(x, ...) standardGeneric("rotate"))}
if (!isGeneric("rescale")) {setGeneric("rescale", function(x, ...) standardGeneric("rescale"))}
if (!isGeneric("select")) {setGeneric("select", function(x, ...) standardGeneric("select"))}
if (!isGeneric("setMinMax")) {setGeneric("setMinMax", function(x, ...) standardGeneric("setMinMax"))}
if (!isGeneric("scale")) {setGeneric("scale", function(x) standardGeneric("scale"))}
if (!isGeneric("shift")) {setGeneric("shift", function(x, ...) standardGeneric("shift"))}
if (!isGeneric("stdev")) { setGeneric("stdev", function(x, ...) standardGeneric("stdev")) }
if (!isGeneric("subset")) {setGeneric("subset", function(x, ...) standardGeneric("subset")) }
if (!isGeneric("summary")) {setGeneric("summary", function(object, ...) standardGeneric("summary")) }
if (!isGeneric("t")) { setGeneric("t", function(x) standardGeneric("t"))}
if (!isGeneric("tail")) { setGeneric("tail", function(x, ...) standardGeneric("tail"))}
if (!isGeneric("text")) { setGeneric("text", function(x, ...) standardGeneric("text")) }
if (!isGeneric("trans")) { setGeneric("trans", function(x, ...) standardGeneric("trans"))}
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
if (!isGeneric("which.lyr")) {setGeneric("which.lyr", function(x) standardGeneric("which.lyr"))}
if (!isGeneric("writeStart")) {	setGeneric("writeStart", function(x, filename, ...)	standardGeneric("writeStart")) }
if (!isGeneric("writeStop")) { setGeneric("writeStop", function(x, ...) standardGeneric("writeStop")) }
if (!isGeneric("writeValues")) { setGeneric("writeValues", function(x, v, ...) standardGeneric("writeValues")) }
if (!isGeneric("writeRaster")) {setGeneric("writeRaster", function(x, filename, ...) standardGeneric("writeRaster"))}
if (!isGeneric("xmin")) {setGeneric("xmin", function(x) standardGeneric("xmin"))}
if (!isGeneric("xmax")) {setGeneric("xmax", function(x)	standardGeneric("xmax"))}
if (!isGeneric("ymin")) {setGeneric("ymin", function(x)	standardGeneric("ymin"))}
if (!isGeneric("ymax")) {setGeneric("ymax", function(x)	standardGeneric("ymax"))}
if (!isGeneric("xmin<-")) { setGeneric("xmin<-", function(x, ..., value) standardGeneric("xmin<-"))}
if (!isGeneric("xmax<-")) { setGeneric("xmax<-", function(x, ..., value) standardGeneric("xmax<-"))}
if (!isGeneric("ymin<-")) { setGeneric("ymin<-", function(x, ..., value) standardGeneric("ymin<-"))}
if (!isGeneric("ymax<-")) { setGeneric("ymax<-", function(x, ..., value) standardGeneric("ymax<-"))}

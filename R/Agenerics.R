#if (!isGeneric("#")) { setGeneric("#", function(x, ...) standardGeneric("#")) }

if (!isGeneric("k_means")) {setGeneric("k_means", function(x, ...) standardGeneric("k_means"))}
if (!isGeneric("princomp")) {setGeneric("princomp", function(x, ...) standardGeneric("princomp"))}
if (!isGeneric("extractRange")) { setGeneric("extractRange", function(x, y, ...) standardGeneric("extractRange"))}

if (!isGeneric("layerCor")) {setGeneric("layerCor", function(x, ...) standardGeneric("layerCor"))}

if (!isGeneric("metags")) {setGeneric("metags", function(x, ...) standardGeneric("metags"))}
if (!isGeneric("metags<-")) {setGeneric("metags<-", function(x, ..., value) standardGeneric("metags<-"))}

if (!isGeneric("forceCCW")) {setGeneric("forceCCW", function(x, ...) standardGeneric("forceCCW"))}
if (!isGeneric("addCats")) {setGeneric("addCats", function(x, ...) standardGeneric("addCats"))}
if (!isGeneric("regress")) {setGeneric("regress", function(y, x, ...) standardGeneric("regress"))}

if (!isGeneric("panel")) {setGeneric("panel", function(x, ...) standardGeneric("panel"))}
#if (!isGeneric("colSums")) {setGeneric("colSums", function(x, ...) standardGeneric("colSums"))}	
#if (!isGeneric("rowSums")) {setGeneric("rowSums", function(x, ...) standardGeneric("rowSums"))}	
#if (!isGeneric("colMeans")) {setGeneric("colMeans", function(x, ...) standardGeneric("colMeans"))}	
#if (!isGeneric("rowMeans")) {setGeneric("rowMeans", function(x, ...) standardGeneric("rowMeans"))}	

if (!isGeneric("logic")) {setGeneric("logic", function(x, ...) standardGeneric("logic"))}
if (!isGeneric("compare")) {setGeneric("compare", function(x, y, ...) standardGeneric("compare"))}

if (!isGeneric("meta")) {setGeneric("meta", function(x, ...) standardGeneric("meta"))}

if (!isGeneric("rangeFill")) {setGeneric("rangeFill", function(x, ...) standardGeneric("rangeFill"))}
if (!isGeneric("roll")) {setGeneric("roll", function(x, ...) standardGeneric("roll"))}
if (!isGeneric("elongate")) {setGeneric("elongate", function(x, ...) standardGeneric("elongate"))}

if (!isGeneric("update")) {setGeneric("update", function(object, ...) standardGeneric("update"))}

if (!isGeneric("viewshed")) {setGeneric("viewshed", function(x, ...) standardGeneric("viewshed"))}
if (!isGeneric("sieve")) {setGeneric("sieve", function(x, ...) standardGeneric("sieve"))}

if (!isGeneric("rasterizeWin")) {setGeneric("rasterizeWin", function(x, y, ...) standardGeneric("rasterizeWin"))}
if (!isGeneric("interpNear")) {setGeneric("interpNear", function(x, y, ...) standardGeneric("interpNear"))}
if (!isGeneric("interpIDW")) {setGeneric("interpIDW", function(x, y, ...) standardGeneric("interpIDW"))}

if (!isGeneric("normalize.longitude")) {setGeneric("normalize.longitude", function(x, ...) standardGeneric("normalize.longitude"))}

if (!isGeneric("allNA")) {setGeneric("allNA", function(x, ...) standardGeneric("allNA"))}
if (!isGeneric("noNA")) {setGeneric("noNA", function(x, ...) standardGeneric("noNA"))}
if (!isGeneric("countNA")) {setGeneric("countNA", function(x, ...) standardGeneric("countNA"))}

if (!isGeneric("scoff")) {setGeneric("scoff", function(x, ...) standardGeneric("scoff"))}
if (!isGeneric("scoff<-")) {setGeneric("scoff<-", function(x, ..., value) standardGeneric("scoff<-"))}

if (!isGeneric("blocks")) {setGeneric("blocks", function(x, ...) standardGeneric("blocks"))}
if (!isGeneric("droplevels")) {setGeneric("droplevels", function(x, ...) standardGeneric("droplevels"))}

if (!isGeneric("str")) { setGeneric("str", function(object, ...) standardGeneric("str"))}

if (!isGeneric("plet")) { setGeneric("plet", function(x, ...) standardGeneric("plet"))}

if (!isGeneric("combineGeoms")) {setGeneric("combineGeoms", function(x, y, ...) standardGeneric("combineGeoms"))}
if (!isGeneric("concats")) {setGeneric("concats", function(x, ...) standardGeneric("concats"))}
if (!isGeneric("has.colors")) {setGeneric("has.colors", function(x, ...) standardGeneric("has.colors"))}
if (!isGeneric("has.RGB")) {setGeneric("has.RGB", function(x, ...) standardGeneric("has.RGB"))}

if (!isGeneric("emptyGeoms")) {setGeneric("emptyGeoms", function(x, ...) standardGeneric("emptyGeoms"))}

if (!isGeneric("serialize")) {setGeneric("serialize", function(object, connection, ascii = FALSE, xdr = TRUE, version = NULL,  refhook = NULL) standardGeneric("serialize"))}

if (!isGeneric("saveRDS")) {setGeneric("saveRDS", function (object, file="", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL) standardGeneric("saveRDS"))}


if (!isGeneric("query")) {setGeneric("query", function(x, ...) standardGeneric("query"))}

if (!isGeneric("set.values")) {setGeneric("set.values", function(x, ...) standardGeneric("set.values"))}
if (!isGeneric("set.ext")) {setGeneric("set.ext", function(x, ...) standardGeneric("set.ext"))}
if (!isGeneric("set.names")) {setGeneric("set.names", function(x, ...) standardGeneric("set.names"))}
if (!isGeneric("set.crs")) {setGeneric("set.crs", function(x, ...) standardGeneric("set.crs"))}
if (!isGeneric("set.RGB")) {setGeneric("set.RGB", function(x, ...) standardGeneric("set.RGB"))}


if (!isGeneric("math")) { setGeneric("math", function(x, ...) standardGeneric("math")) }

if (!isGeneric("all.equal")) { setGeneric("all.equal", function(target, current, ...) standardGeneric("all.equal")) }
if (!isGeneric("impose")) { setGeneric("impose", function(x, ...) standardGeneric("impose")) }
if (!isGeneric("density")) { setGeneric("density", function(x, ...) standardGeneric("density"))}
if (!isGeneric("densify")) { setGeneric("densify", function(x, ...) standardGeneric("densify"))}


if (!isGeneric("selectHighest")) {setGeneric("selectHighest", function(x, ...) standardGeneric("selectHighest"))}
if (!isGeneric("focal3D")) { setGeneric("focal3D", function(x, ...) standardGeneric("focal3D")) }
if (!isGeneric("focalReg")) { setGeneric("focalReg", function(x, ...) standardGeneric("focalReg")) }
if (!isGeneric("focalCpp")) { setGeneric("focalCpp", function(x, ...) standardGeneric("focalCpp")) }
if (!isGeneric("focalPairs")) { setGeneric("focalPairs", function(x, ...) standardGeneric("focalPairs")) }
if (!isGeneric("focalCor")) { setGeneric("focalCor", function(x, ...) standardGeneric("focalCor")) }

if (!isGeneric("clearance")) {setGeneric("clearance", function(x, ...) standardGeneric("clearance"))}
if (!isGeneric("width")) {setGeneric("width", function(x, ...) standardGeneric("width"))}
if (!isGeneric("simplifyGeom")) {setGeneric("simplifyGeom", function(x, ...) standardGeneric("simplifyGeom"))}
if (!isGeneric("thinGeom")) {setGeneric("thinGeom", function(x, ...) standardGeneric("thinGeom"))}
if (!isGeneric("mergeLines")) {setGeneric("mergeLines", function(x, ...) standardGeneric("mergeLines"))}
if (!isGeneric("mergeTime")) {setGeneric("mergeTime", function(x, ...) standardGeneric("mergeTime"))}
if (!isGeneric("fillTime")) {setGeneric("fillTime", function(x, ...) standardGeneric("fillTime"))}
if (!isGeneric("makeNodes")) {setGeneric("makeNodes", function(x, ...) standardGeneric("makeNodes"))}
if (!isGeneric("removeDupNodes")) {setGeneric("removeDupNodes", function(x, ...) standardGeneric("removeDupNodes"))}
if (!isGeneric("snap")) {setGeneric("snap", function(x, ...) standardGeneric("snap"))}

if (!isGeneric("weighted.mean")) {setGeneric("weighted.mean", function(x, w, ...) standardGeneric("weighted.mean"))}

if (!isGeneric("split")) {setGeneric("split", function(x, f, drop = FALSE, ...) standardGeneric("split"))}
if (!isGeneric("cellSize")) {setGeneric("cellSize", function(x, ...) standardGeneric("cellSize"))}

if (!isGeneric("na.omit")) {setGeneric("na.omit", function(object, ...) standardGeneric("na.omit"))}
if (!isGeneric("catalyze")) {setGeneric("catalyze", function(x, ...) standardGeneric("catalyze"))}
if (!isGeneric("activeCat")) {setGeneric("activeCat", function(x, ...) standardGeneric("activeCat"))}
if (!isGeneric("activeCat<-")) {setGeneric("activeCat<-", function(x, ..., value) standardGeneric("activeCat<-"))}

if (!isGeneric("sharedPaths")) {setGeneric("sharedPaths", function(x, ...) standardGeneric("sharedPaths"))}
if (!isGeneric("isTRUE")) {setGeneric("isTRUE", function(x) standardGeneric("isTRUE"))}
if (!isGeneric("subst")) {setGeneric("subst", function(x, ...) standardGeneric("subst"))}

if (!isGeneric("colorize")) {setGeneric("colorize", function(x, ...) standardGeneric("colorize"))}

if (!isGeneric("RGB")) {setGeneric("RGB", function(x, ...) standardGeneric("RGB"))}
if (!isGeneric("RGB<-")) {setGeneric("RGB<-", function(x, ..., value) standardGeneric("RGB<-"))}

if (!isGeneric("autocor")) {setGeneric("autocor", function(x, ...) standardGeneric("autocor"))}
if (!isGeneric("delaunay")) {setGeneric("delaunay", function(x, ...) standardGeneric("delaunay"))}
if (!isGeneric("voronoi")) {setGeneric("voronoi", function(x, ...) standardGeneric("voronoi"))}
if (!isGeneric("convHull")) {setGeneric("convHull", function(x, ...) standardGeneric("convHull"))}
if (!isGeneric("minRect")) {setGeneric("minRect", function(x, ...) standardGeneric("minRect"))}
if (!isGeneric("minCircle")) {setGeneric("minCircle", function(x, ...) standardGeneric("minCircle"))}

#if (!isGeneric("which.related")) {setGeneric("which.related", function(x, y, ...) standardGeneric("which.related"))}
if (!isGeneric("is.related")) {setGeneric("is.related", function(x, y, ...) standardGeneric("is.related"))}
if (!isGeneric("relate")) {setGeneric("relate", function(x, y, ...) standardGeneric("relate"))}
if (!isGeneric("intersect")) {setGeneric("intersect", function(x, y) standardGeneric("intersect"))}

if (!isGeneric("not.na")) {setGeneric("not.na", function(x, ...) standardGeneric("not.na"))}

if (!isGeneric("erase")) {setGeneric("erase", function(x, y, ...) standardGeneric("erase"))}
if (!isGeneric("gaps")) {setGeneric("gaps", function(x, ...) standardGeneric("gaps"))}
if (!isGeneric("is.rotated")) {setGeneric("is.rotated", function(x, ...) standardGeneric("is.rotated"))}

if (!isGeneric("is.int")) {setGeneric("is.int", function(x) standardGeneric("is.int"))}
if (!isGeneric("as.int")) {setGeneric("as.int", function(x, ...) standardGeneric("as.int"))}
if (!isGeneric("is.bool")) {setGeneric("is.bool", function(x) standardGeneric("is.bool"))}
if (!isGeneric("as.bool")) {setGeneric("as.bool", function(x, ...) standardGeneric("as.bool"))}
if (!isGeneric("nearby")) {setGeneric("nearby", function(x, ...) standardGeneric("nearby"))}
if (!isGeneric("nearest")) {setGeneric("nearest", function(x, ...) standardGeneric("nearest"))}
if (!isGeneric("cartogram")) {setGeneric("cartogram", function(x, ...) standardGeneric("cartogram"))}
if (!isGeneric("dots")) {setGeneric("dots", function(x, ...) standardGeneric("dots"))}
if (!isGeneric("crds")) {setGeneric("crds", function(x, ...) standardGeneric("crds"))}
if (!isGeneric("symdif")) {setGeneric("symdif", function(x, y, ...) standardGeneric("symdif"))}
if (!isGeneric("median")) {setGeneric("median", function(x, na.rm) standardGeneric("median"))}
if (!isGeneric("polys")) {setGeneric("polys", function(x,...) standardGeneric("polys"))}
if (!isGeneric("centroids")) {setGeneric("centroids", function(x, ...) standardGeneric("centroids"))}
if (!isGeneric("coltab")) {setGeneric("coltab", function(x, ...) standardGeneric("coltab"))}
if (!isGeneric("coltab<-")) { setGeneric("coltab<-", function(x, ..., value) standardGeneric("coltab<-")) }
if (!isGeneric("deepcopy")) { setGeneric("deepcopy", function(x, ...) standardGeneric("deepcopy")) }
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
if (!isGeneric("makeValid")) {setGeneric("makeValid", function(x,...) standardGeneric("makeValid"))}
if (!isGeneric("is.valid")) {setGeneric("is.valid", function(x,...) standardGeneric("is.valid"))}
if (!isGeneric("is.empty")) {setGeneric("is.empty", function(x,...) standardGeneric("is.empty"))}
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
if (!isGeneric("hasMinMax")) {setGeneric("hasMinMax", function(x) standardGeneric("hasMinMax"))}
if (!isGeneric("minmax")) {setGeneric("minmax", function(x, ...) standardGeneric("minmax"))}
if (!isGeneric("nsrc")) { setGeneric("nsrc", function(x) standardGeneric("nsrc")) }
if (!isGeneric("perim")) {setGeneric("perim", function(x, ...) standardGeneric("perim"))}
if (!isGeneric("project")) {setGeneric("project", function(x,...) standardGeneric("project"))}
if (!isGeneric("wrapCache")) {setGeneric("wrapCache", function(x, ...) standardGeneric("wrapCache"))}
if (!isGeneric("wrap")) {setGeneric("wrap", function(x, ...) standardGeneric("wrap"))}
if (!isGeneric("unwrap")) {setGeneric("unwrap", function(x, ...) standardGeneric("unwrap"))}
if (!isGeneric("cats")) { setGeneric("cats", function(x, ...) standardGeneric("cats")) }
if (!isGeneric("categories")) { setGeneric("categories", function(x, ...) standardGeneric("categories")) }
if (!isGeneric("set.cats")) { setGeneric("set.cats", function(x, ...) standardGeneric("set.cats")) }
if (!isGeneric("as.raster")) { setGeneric("as.raster", function(x, ...) standardGeneric("as.raster"))}
if (!isGeneric("rast") ) { setGeneric("rast", function(x, ...) standardGeneric("rast")) }
if (!isGeneric("rev") ) { setGeneric("rev", function(x) standardGeneric("rev")) }
if (!isGeneric("sds") ) { setGeneric("sds", function(x, ...) standardGeneric("sds")) }
if (!isGeneric("sprc") ) { setGeneric("sprc", function(x, ...) standardGeneric("sprc")) }
if (!isGeneric("svc") ) { setGeneric("svc", function(x, ...) standardGeneric("svc")) }
if (!isGeneric("sel")) {setGeneric("sel", function(x, ...) standardGeneric("sel"))}
if (!isGeneric("segregate")) {setGeneric("segregate", function(x, ...) standardGeneric("segregate"))}
if (!isGeneric("selectRange")) {setGeneric("selectRange", function(x, ...) standardGeneric("selectRange"))}
if (!isGeneric("setValues")) {setGeneric("setValues", function(x, values, ...) standardGeneric("setValues"))}
if (!isGeneric("expanse")) {setGeneric("expanse", function(x, ...) standardGeneric("expanse"))}
if (!isGeneric("size")) {setGeneric("size", function(x, ...) standardGeneric("size"))}
if (!isGeneric("inMemory")) {setGeneric("inMemory", function(x, ...) standardGeneric("inMemory"))}
if (!isGeneric("sources")) {setGeneric("sources", function(x, ...) standardGeneric("sources"))}
if (!isGeneric("spatSample")) { setGeneric("spatSample", function(x, ...) standardGeneric("spatSample"))}
if (!isGeneric("terrain")) {setGeneric("terrain", function(x, ...) standardGeneric("terrain"))}
if (!isGeneric("has.time")) {setGeneric("has.time", function(x,...) standardGeneric("has.time"))}
if (!isGeneric("time")) {setGeneric("time", function(x,...) standardGeneric("time"))}
if (!isGeneric("time<-")) {setGeneric("time<-", function(x, ..., value) standardGeneric("time<-"))}
if (!isGeneric("timeInfo")) {setGeneric("timeInfo", function(x,...) standardGeneric("timeInfo"))}
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
if (!isGeneric("approximate")) {setGeneric("approximate", function(x, ...) standardGeneric("approximate"))}

if (!isGeneric("as.data.frame")) { setGeneric("as.data.frame", function(x, row.names = NULL, optional = FALSE, ...) standardGeneric("as.data.frame")) }

if (!isGeneric("as.list")) { setGeneric("as.list", function(x, ...) standardGeneric("as.list"))}
if (!isGeneric("as.factor")) {setGeneric("as.factor", function(x) standardGeneric("as.factor"))}

if (!isGeneric("atan2")) { setGeneric("atan2", function(y, x) standardGeneric("atan2"))}
if (!isGeneric("atan_2")) { setGeneric("atan_2", function(y, x, ...) standardGeneric("atan_2"))}
if (!isGeneric("barplot")) {setGeneric("barplot", function(height,...) standardGeneric("barplot"))}
#if (!isGeneric("bndbox")) {setGeneric("bndbox", function(obj) standardGeneric("bndbox"))}
if (!isGeneric("boundaries")) {	setGeneric("boundaries", function(x, ...) standardGeneric("boundaries"))}
if (!isGeneric("boxplot")) { setGeneric("boxplot", function(x, ...) standardGeneric("boxplot"))}
if (!isGeneric("buffer")) {setGeneric("buffer", function(x, ...) standardGeneric("buffer"))}
if (!isGeneric("clamp")) { setGeneric("clamp", function(x, ...) standardGeneric("clamp")) }
if (!isGeneric("clamp_ts")) { setGeneric("clamp_ts", function(x, ...) standardGeneric("clamp_ts")) }
if (!isGeneric("click")) {setGeneric("click", function(x, ...)standardGeneric("click"))}
if (!isGeneric("contour")) { setGeneric("contour", function(x,...) standardGeneric("contour"))}
if (!isGeneric("cover")) {setGeneric("cover", function(x, y, ...) standardGeneric("cover"))}
if (!isGeneric("crop")) { setGeneric("crop", function(x, y, ...) standardGeneric("crop")) }
if (!isGeneric("crs")) { setGeneric("crs", function(x, ...)	standardGeneric("crs")) }
if (!isGeneric("crs<-")) { setGeneric("crs<-", function(x, ..., value) standardGeneric("crs<-")) }
if (!isGeneric("density")) { setGeneric("density", function(x, ...) standardGeneric("density"))}
if (!isGeneric("aggregate")) {setGeneric("aggregate", function(x, ...) standardGeneric("aggregate"))}
if (!isGeneric("disagg")) {setGeneric("disagg", function(x, ...) standardGeneric("disagg"))}
#if (!isGeneric("costDistance")) {setGeneric("costDistance", function(x, ...)standardGeneric("costDistance"))}
if (!isGeneric("gridDistance")) {setGeneric("gridDistance", function(x, ...)standardGeneric("gridDistance"))}
if (!isGeneric("costDist")) {setGeneric("costDist", function(x, ...)standardGeneric("costDist"))}
if (!isGeneric("gridDist")) {setGeneric("gridDist", function(x, ...)standardGeneric("gridDist"))}
if (!isGeneric("distance")) {setGeneric("distance", function(x, y, ...)standardGeneric("distance"))}
if (!isGeneric("direction")) {setGeneric("direction", function(x, ...)standardGeneric("direction"))}
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
if (!isGeneric("inext")) {setGeneric("inext", function(x, ...) standardGeneric("inext"))}
if (!isGeneric("interpolate")) { setGeneric("interpolate", function(object, ...) standardGeneric("interpolate"))}
if (!isGeneric("is.factor")) {setGeneric("is.factor", function(x) standardGeneric("is.factor"))}
if (!isGeneric("is.lonlat")) { setGeneric("is.lonlat", function(x, ...) standardGeneric("is.lonlat"))}
if (!isGeneric("mask")) { setGeneric("mask", function(x, mask, ...) standardGeneric("mask")) }
if (!isGeneric("match")) { setGeneric("match", function(x, table, nomatch=NA_integer_, incomparables=NULL)		standardGeneric("match"))}
if (!isGeneric("modal")) {setGeneric("modal", function(x, ...) standardGeneric("modal"))}
if (!isGeneric("mosaic")) {setGeneric("mosaic", function(x, y, ...) standardGeneric("mosaic"))}
if (!isGeneric("ncell")) { setGeneric("ncell", function(x) standardGeneric("ncell")) }
if (!isGeneric("nrow")) { setGeneric("nrow", function(x) standardGeneric("nrow")) }
if (!isGeneric("ncol")) { setGeneric("nrow", function(x) standardGeneric("nrow")) }
if (!isGeneric("ncol<-")) { setGeneric("ncol<-", function(x, ..., value) standardGeneric("ncol<-")) }
if (!isGeneric("nrow<-")) { setGeneric("nrow<-", function(x, ..., value) standardGeneric("nrow<-")) }
if (!isGeneric("origin")) {	setGeneric("origin", function(x, ...) standardGeneric("origin")) }
if (!isGeneric("origin<-")) {setGeneric("origin<-", function(x, value)	standardGeneric("origin<-"))}
if (!isGeneric("pairs")) { setGeneric("pairs", function(x, ...)	standardGeneric("pairs"))}
if (!isGeneric("patches")) {setGeneric("patches", function(x, ...) standardGeneric("patches"))}
if (!isGeneric("persp")) { setGeneric("persp", function(x,...) standardGeneric("persp")) }
if (!isGeneric("plot")) { setGeneric("plot", function(x, y,...) standardGeneric("plot"))}
if (!isGeneric("plotRGB")) { setGeneric("plotRGB", function(x, ...)standardGeneric("plotRGB"))}
if (!isGeneric("predict")) {setGeneric("predict", function(object, ...) standardGeneric("predict"))}
if (!isGeneric("quantile")) {setGeneric("quantile", function(x, ...)standardGeneric("quantile"))}
if (!isGeneric("rasterize")) {setGeneric("rasterize", function(x, y, ...) standardGeneric("rasterize"))}
if (!isGeneric("rasterizeGeom")) {setGeneric("rasterizeGeom", function(x, y, ...) standardGeneric("rasterizeGeom"))}
if (!isGeneric("readStart")) {setGeneric("readStart", function(x, ...) standardGeneric("readStart"))}
if (!isGeneric("readStop")) {setGeneric("readStop", function(x)	standardGeneric("readStop"))}
if (!isGeneric("res")) { setGeneric("res", function(x) standardGeneric("res")) }
if (!isGeneric("res<-")) { setGeneric("res<-", function(x, value) standardGeneric("res<-")) }
if (!isGeneric("rectify")) {setGeneric("rectify", function(x, ...) standardGeneric("rectify"))}
if (!isGeneric("resample")) { setGeneric("resample", function(x, y, ...) standardGeneric("resample"))}
if (!isGeneric("spin")) {setGeneric("spin", function(x, ...) standardGeneric("spin"))}
if (!isGeneric("rotate")) {setGeneric("rotate", function(x, ...) standardGeneric("rotate"))}
if (!isGeneric("rescale")) {setGeneric("rescale", function(x, ...) standardGeneric("rescale"))}
#if (!isGeneric("select")) {setGeneric("select", function(x, ...) standardGeneric("select"))}
if (!isGeneric("setMinMax")) {setGeneric("setMinMax", function(x, ...) standardGeneric("setMinMax"))}
if (!isGeneric("scale")) {setGeneric("scale", function(x, center=TRUE, scale=TRUE) standardGeneric("scale"))}
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
if (!isGeneric("rcl")) { setGeneric("rcl", function(x, ...) standardGeneric("rcl")) }

if (!isGeneric("yFromRow")) { setGeneric("yFromRow", function(object, row) standardGeneric("yFromRow")) }
if (!isGeneric("xFromCol")) { setGeneric("xFromCol", function(object, col) standardGeneric("xFromCol")) }
if (!isGeneric("colFromX")) { setGeneric("colFromX", function(object, x) standardGeneric("colFromX")) }
if (!isGeneric("rowFromY")) { setGeneric("rowFromY", function(object, y) standardGeneric("rowFromY")) }
if (!isGeneric("cellFromXY")) { setGeneric("cellFromXY", function(object, xy) standardGeneric("cellFromXY")) }
if (!isGeneric("cellFromRowCol")) { setGeneric("cellFromRowCol", function(object, row, col, ...) standardGeneric("cellFromRowCol")) }
if (!isGeneric("rowColCombine")) { setGeneric("rowColCombine", function(object, row, col, ...) standardGeneric("rowColCombine")) }
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
if (!isGeneric("union")) {setGeneric("union", function(x, y)standardGeneric("union"))}
if (!isGeneric("unique")) { setGeneric("unique", function(x, incomparables=FALSE, ...) standardGeneric("unique")) }
if (!isGeneric("values")) { setGeneric("values", function(x, ...) standardGeneric("values")) }
if (!isGeneric("values<-")) { setGeneric("values<-", function(x, value) standardGeneric("values<-"))}
if (!isGeneric("where.max")) {setGeneric("where.max", function(x, ...) standardGeneric("where.max"))}
if (!isGeneric("where.min")) {setGeneric("where.min", function(x, ...) standardGeneric("where.min"))}
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
if (!isGeneric("zoom")) {setGeneric("zoom", function(x, ...)standardGeneric("zoom"))}



setClass("PackedSpatVector",
	representation (
		type = "character",
		crs = "character",
		coordinates = "matrix",
		index = "matrix",
		attributes = "data.frame"
	),
	prototype (
		type= "",
		crs = ""
	)
)


setClass("PackedSpatRaster",
	representation (
		definition = "character",
		values = "matrix",
		attributes = "list"
	),
	prototype (
		attributes = list()
	)
)


.packVector <- function(x) {
	vd <- methods::new("PackedSpatVector")
	vd@type <- geomtype(x)
	vd@crs <- as.character(crs(x))
	stopifnot(vd@type %in% c("points", "lines", "polygons"))
	g <- geom(x)
	vd@coordinates <- g[, c("x", "y")]
	j <- c(1,2, grep("hole", colnames(g)))
	g <- g[,j]
	i <- which(!duplicated(g))
	vd@index <- cbind(g[i, ], start=i)
	vd
}

setMethod("wrap", signature(x="Spatial"), 
	function(x) {
		pv <- .packVector(x)
		if (methods::.hasSlot(x, "data")) {
			pv@attributes <- x@data
		}
		pv
	}
)


setMethod("wrap", signature(x="SpatVector"), 
	function(x) {
		pv <- .packVector(x)
		pv@attributes <- as.data.frame(x)
		pv
	}
)


setMethod("vect", signature(x="PackedSpatVector"), 
	function(x) {
		p <- methods::new("SpatVector")
		p@ptr <- SpatVector$new()
		if (!is.na(x@crs)) {
			crs(p) <- x@crs
		}
		if (nrow(x@coordinates) == 0) {
			return(p)
		}

		n <- ncol(x@index)
		reps <- diff(c(x@index[,n], nrow(x@coordinates)+1))
		i <- rep(1:nrow(x@index), reps)
		if (n == 2) { 
			p@ptr$setGeometry(x@type, x@index[i,1], x@index[i,2], x@coordinates[,1], x@coordinates[,2], rep(0, nrow(x@coordinates)))
		} else {
			p@ptr$setGeometry(x@type, x@index[i,1], x@index[i,2], x@coordinates[,1], x@coordinates[,2], x@index[i,3])
		} 
		if (nrow(x@attributes) > 0) {
			values(p) <- x@attributes
		}
		messages(p, "pack")
	}
)

setMethod("show", signature(object="PackedSpatVector"), 
	function(object) {
		print(paste("This is a", class(object), "object. Use 'terra::vect()' to unpack it"))
	}
)



setMethod("as.character", signature(x="SpatRaster"), 
	function(x) {
		e <- as.vector(ext(x))
		d <- crs(x, describe=TRUE)
		if (!(is.na(d$authority) || is.na(d$code))) {
			crs <- paste0(", crs='", d$authority, ":", d$code, "'")
		} else {
			d <- crs(x)
			crs <- ifelse(d=="", ", crs=''", paste0(", crs='", d, "'"))
			crs <- gsub("\n[ ]+", "", crs)
		}
		nms <- paste0(", names=c('", paste(names(x), collapse="', '"), "')")
		paste0("rast(", 
				"ncols=", ncol(x),
				", nrows=", nrow(x),
				", nlyrs=", nlyr(x),
				", xmin=",e[1],
				", xmax=",e[2],
				", ymin=",e[3],
				", ymax=",e[4],
				nms, 
				crs, ")" 
		)
	}
)
#eval(parse(text=as.character(s)))


setMethod("wrap", signature(x="SpatRaster"), 
	function(x) {
		r <- methods::new("PackedSpatRaster")
		r@definition <- as.character(x)
		r@values <- values(x)
		if (any(is.factor(x))) {
			r@attributes$levels <- levels(x)
		} 
		v <- time(x)
		if (any(!is.na(v))) {
			r@attributes$time <- v
		} 
		v <- units(x)
		if (any(!is.na(v))) {
			r@attributes$units <- v
		} 
		v <- depth(x)
		if (any(!is.na(v))) {
			r@attributes$depth <- v
		} 
		r
	}
)


setMethod("rast", signature(x="PackedSpatRaster"), 
	function(x) {
		r <- eval(parse(text=x@definition))
		values(r) <- x@values
		if (length(x@attributes) > 0) {
			nms <- names(x@attributes)
			if (all(nms %in% c("levels", "time", "units", "depth"))) {
				time(r) <- x@attributes$time
				units(r) <- x@attributes$units
				levels(r) <- x@attributes$levels
				depth(r) <- x@attributes$depth
			} else {
				levels(r) <- x@attributes
			}
		}
		r
	}
)

setMethod("show", signature(object="PackedSpatRaster"), 
	function(object) {
		print(paste("This is a", class(object), "object. Use 'terra::rast()' to unpack it"))
	}
)



setMethod("serialize", signature(object="SpatVector"), 
	function(object, connection, ascii = FALSE, xdr = TRUE, version = NULL, refhook = NULL) {
		object = wrap(object)
		serialize(object, connection=connection, ascii = ascii, xdr = xdr, version = version, refhook = refhook)
	}
)


setMethod("saveRDS", signature(object="SpatVector"), 
	function(object, file="", ascii = FALSE, version = NULL, compress=TRUE, refhook = NULL) {
		object = wrap(object)
		saveRDS(object, file=file, ascii = ascii, version = version, compress=compress, refhook = refhook)
	}
)


setMethod("serialize", signature(object="SpatRaster"), 
	function(object, connection, ascii = FALSE, xdr = TRUE, version = NULL, refhook = NULL) {
		if (!all(inMemory(object))) {
			opt <- spatOptions()
			if (object@ptr$canProcessInMemory(opt)) {
				set.values(object)
			} else {
				error("Cannot be loaded into memory which is required for serialize")
			}
		}
		object <- wrap(object)
		serialize(object, connection=connection, ascii = ascii, xdr = xdr, version = version, refhook = refhook)
	}
)


setMethod("saveRDS", signature(object="SpatRaster"), 
	function(object, file="", ascii = FALSE, version = NULL, compress=TRUE, refhook = NULL) {
		if (!all(inMemory(object))) {
			opt <- spatOptions()
			if (object@ptr$canProcessInMemory(opt)) {
				set.values(object)
			} else {
				error("Cannot be loaded into memory which is required for saveRDS")
			}
		}
		object = wrap(object)
		saveRDS(object, file=file, ascii = ascii, version = version, compress=compress, refhook = refhook)
	}
)

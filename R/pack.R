

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
		attributes = "data.frame"
	),	
)


.packVector <- function(x) {
	vd <- methods::new("PackedSpatVector")
	vd@type <- geomtype(x)
	vd@crs <- as.character(crs(x))
	stopifnot(vd@type %in% c("points", "lines", "polygons"))
	g <- as.matrix(geom(x))
	vd@coordinates <- g[, c("x", "y")]
	j <- c(1,2, grep("hole", colnames(g)))
	g <- g[,j]
	i <- which(!duplicated(g))
	vd@index <- cbind(g[i, ], start=i)
	vd
}

setMethod("pack", signature(x="Spatial"), 
	function(x, ...) {
		pv <- .packVector(x)
		if (methods::.hasSlot(x, "data")) {
			pv@attributes <- x@data	
		}
		pv
	}
)


setMethod("pack", signature(x="SpatVector"), 
	function(x, ...) {
		pv <- .packVector(x)
		pv@attributes <- as.data.frame(x)
		pv
	}
)


setMethod("vect", signature(x="PackedSpatVector"), 
	function(x, ...) {
		p <- methods::new("SpatVector")
		p@ptr <- SpatVector$new()
		if (!is.na(x@crs)) {
			p@ptr$crs <- x@crs
		}
		if (nrow(x@coordinates) == 0) {
			return(p)
		}

		n <- ncol(x@index)
		reps <- diff(c(x@index[,n], nrow(x@coordinates)+1))
		i <- rep(1:nrow(x@index), reps)
		if (n == 3) { 
			p@ptr$setGeometry(x@type, x@index[i,1], x@index[i,2], x@coordinates[,1], x@coordinates[,2], rep(0, nrow(x@coordinates)))		
		} else {
			p@ptr$setGeometry(x@type, x@index[i,1], x@index[i,2], x@coordinates[,1], x@coordinates[,2], x@index[i,3])
		} 
		if (nrow(x@attributes) > 0) {
			values(p) <- x@attributes
		}
		show_messages(p)
	}
)

setMethod("show", signature(object="PackedSpatVector"), 
	function(object) {
		print(paste("This is a", class(object), "object. Use 'vect()' to unpack it"))
	}
)





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
		show_messages(p, "pack")
	}
)

setMethod("show", signature(object="PackedSpatVector"), 
	function(object) {
		print(paste("This is a", class(object), "object. Use 'vect()' to unpack it"))
	}
)





setMethod("as.character", signature(x="SpatRaster"), 
	function(x, ...) {
		e <- as.vector(ext(x))
		crs <- crs(x)
		crs <- ifelse(is.na(crs), ", crs=''", paste0(", crs='", crs, "'"))
		crs <- gsub("\n[ ]+", "", crs)
		paste0("rast(", 
				"ncol=", ncol(x),
				", nrow=", nrow(x),
				", nlyr=", nlyr(x),
				", xmin=",e[1],
				", xmax=",e[2],
				", ymin=",e[3],
				", ymax=",e[4],
				crs, ")" 
		)
	}
)
#eval(parse(text=as.character(raster())))
#eval(parse(text=as.character(stack())))


setMethod("pack", signature(x="SpatRaster"), 
	function(x, ...) {
		r <- methods::new("PackedSpatRaster")
		r@definition = as.character(x)
		r@values = values(x)
		#r@attributes = "data.frame"
		r
	}
)


setMethod("rast", signature(x="PackedSpatRaster"), 
	function(x, ...) {
		r <- eval(parse(text=x@definition))
		values(r) <- x@values
		r
	}
)

setMethod("show", signature(object="PackedSpatRaster"), 
	function(object) {
		print(paste("This is a", class(object), "object. Use 'rast()' to unpack it"))
	}
)

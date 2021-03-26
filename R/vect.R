
#setMethod("row.names", signature(x="SpatVector"), 
#	function(x) {
#		1:nrow(x)
#	}
#)

setMethod("vect", signature(x="missing"), 
	function(...) {
		p <- methods::new("SpatVector")
		p@ptr <- SpatVector$new()
		messages(p, "vect")
		return(p)
	}
)

setMethod("vect", signature(x="character"), 
	function(x, ...) {
		p <- methods::new("SpatVector")
		s <- substr(x[1], 1, 5)
		if (s %in% c("POINT", "MULTI", "LINES", "POLYG")) {
#		if (all(grepl("\\(", x) & grepl("\\)", x))) {
			x <- gsub("\n", "", x)
			p@ptr <- SpatVector$new(x)
			dots <- list(...)
			if (!is.null(dots$crs)) {
				crs(p) <- dots$crs
			}
		} else {
			p@ptr <- SpatVector$new()
			x <- normalizePath(x)
			p@ptr$read(x)
		}
		messages(p, "vect")
	}
)


setMethod("vect", signature(x="Spatial"), 
	function(x, ...) {
		methods::as(x, "SpatVector")
	}
)

setMethod("vect", signature(x="sf"), 
	function(x) {
		methods::as(x, "SpatVector")
	}
)

setMethod("vect", signature(x="sfc"), 
	function(x) {
		methods::as(x, "SpatVector")
	}
)

setMethod("vect", signature(x="XY"), #sfg
	function(x) {
		v <- try(.from_sfg(x), silent=TRUE)
		if (inherits(v, "try-error")) {
			error("as,sfg", "coercion failed. You can try coercing via a Spatial* (sp) class")
		} 
		v
	}
)



.checkXYnames <- function(x, warn=FALSE) {
	if (is.null(x)) return(TRUE)
	if (length(x) != 2) {
		error("vect", "coordinate matrix should have 2 columns")
	}
	z <- tolower(x[1:2])
	x <- substr(z, 1, 3)
	y <- substr(x, 1, 1)
	if ((y[1] == "x") & (y[2] == "y")) return(TRUE)
	if ((x[1] == "eas") & (x[2] == "nor")) return(TRUE)
	if ((x[1] == "lon") & (x[2] == "lat")) return(TRUE)
	if (grepl("lon", z[1]) & grepl("lat", z[2])) return(TRUE)

	if ((x[1] == "lat") | (x[2] == "lon")) {
		error("vect", "longitude/latitude in the wrong order")
	} else if ((y[1] == "y") | (y[2] == "x")) {
		error("vect", "x/y in the wrong order")
	} else if ((x[1] == "nor") | (x[2] == "eas")) {
		error("vect", "easting/northing in the wrong order")
	} else if (warn) {
		warn("coordinate names not recognized. Expecting lon/lat, x/y, or easting/northing")
	}
}

setMethod("vect", signature(x="matrix"), 
	function(x, type="points", atts=NULL, crs="", ...) {
		type <- tolower(type)
		stopifnot(type %in% c("points", "lines", "polygons"))
		stopifnot(NCOL(x) > 1)

		p <- methods::new("SpatVector")
		p@ptr <- SpatVector$new()
		nr <- nrow(x)
		if (nr == 0) {
			return(p)
		}

		if (ncol(x) == 2) { 
			.checkXYnames(colnames(x))
			if (type == "points") {	# treat as unique points
				p@ptr$setGeometry(type, 1:nr, rep(1, nr), x[,1], x[,2], rep(FALSE, nr))
			} else {
				p@ptr$setGeometry(type, rep(1, nr), rep(1, nr), x[,1], x[,2], rep(FALSE, nr))
			}
		} else if (ncol(x) == 4) {
			#.checkXYnames(colnames(x)[3:4])
			p@ptr$setGeometry(type, x[,1], x[,2], x[,3], x[,4], rep(FALSE, nr))
		} else if (ncol(x) == 5) {
			#.checkXYnames(colnames(x)[3:4])
			p@ptr$setGeometry(type, x[,1], x[,2], x[,3], x[,4], x[,5])
		} else {
			error("vect", "not an appropriate matrix")
		}
		if (!is.null(atts)) {
			if ((nrow(atts) == nrow(p)) & (ncol(atts) > 0)) {
				values(p) <- atts
			}
		}
		crs(p) <- ifelse(is.na(crs), "", as.character(crs))
		messages(p, "vect")
	}
)


setMethod("$", "SpatVector",  function(x, name) { 
	s <- .subset_cols(x, name, drop=TRUE) 
	s[,1,drop=TRUE]
})


setMethod("[[", c("SpatVector", "numeric", "missing"),
function(x, i, j, ... ,drop=FALSE) {
	s <- .subset_cols(x, i, ..., drop=TRUE)
	s[,,drop=drop]
})


setMethod("[[", c("SpatVector", "character", "missing"),
function(x, i, j, ... ,drop=FALSE) {
	s <- .subset_cols(x, i, ..., drop=TRUE)
	s[,,drop=drop]
})



setReplaceMethod("[[", c("SpatVector", "character", "missing"),
	function(x, i, j, value) {
		
		if (is.null(value)) {
			for (name in i) {
				if (name %in% names(x)) {
					x@ptr$remove_column(name)
				}
			}
			return(x);
		}
		if (length(i) > 1) {
			error("[[<-", "you can only set one variable at a time")
		}
	
		name <- i[1]
		value <- rep(value, length.out=nrow(x))

		if (name %in% names(x)) {
			d <- values(x)
			d[[name]] <- value
			values(x) <- d
		} else {
			if (is.integer(value)) {
				ok <- x@ptr$add_column_long(value, name)
			} else if (is.numeric(value)) {
				ok <- x@ptr$add_column_double(value, name)
			} else {
				ok <- x@ptr$add_column_string(as.character(value), name)
			}
			if (!ok) {
				error("[[<-,SpatVector", "cannot set these values")
			}
		} 
		x
	}
)

setReplaceMethod("[[", c("SpatVector", "numeric", "missing"),
	function(x, i, j, value) {
		stopifnot(i > 0 && i <= ncol(x))
		vn <- names(x)[i]
		x[[vn]] <- value
		x
	}
)



setMethod("$<-", "SpatVector",  
	function(x, name, value) {
		x[[name]] <- value
		x
	}
)



setMethod("vect", signature(x="data.frame"), 
	function(x, geom=c("lon", "lat"), crs=NA, ...) {
		if (length(geom) == 2) {
			v <- vect(as.matrix(x[,geom]), crs=crs)
		} else if (length(geom) == 1) {
			v <- vect(unlist(x[,geom]), crs=crs)
		}
		values(v) <- x
		v
	}
)

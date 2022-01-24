
#setMethod("row.names", signature(x="SpatVector"), 
#	function(x) {
#		1:nrow(x)
#	}
#)

character_crs <- function(crs, caller="") {
	if (is.na(crs)) {
		""
	} else if (!inherits(crs, "character")) {
		warn(caller, "argument 'crs' should be a character value")
		as.character(crs)
	} else {
		crs
	}
}


setMethod("as.vector", signature(x="SpatVector"), 
	function(x, mode="any") {
		if (nrow(x) > 0) {
			lapply(1:nrow(x), function(i) x[i,])
		} else {
			x
		}
	}
)

setMethod("vect", signature(x="missing"), 
	function(x) {
		p <- methods::new("SpatVector")
		p@ptr <- SpatVector$new()
		messages(p, "vect")
		return(p)
	}
)

setMethod("vect", signature(x="character"), 
	function(x, layer="", query="", extent=NULL, filter=NULL, crs="", proxy=FALSE) {
		p <- methods::new("SpatVector")
		s <- substr(x[1], 1, 5)
		if (s %in% c("POINT", "MULTI", "LINES", "POLYG")) {
#		if (all(grepl("\\(", x) & grepl("\\)", x))) {
			x <- gsub("\n", "", x)
			p@ptr <- SpatVector$new(x)
			crs(p) <- crs
		} else {
			p@ptr <- SpatVector$new()
			x <- normalizePath(x)
			proxy <- isTRUE(proxy)
			if (proxy) query <- ""
			if (proxy || is.null(filter)) {
				filter <- vect()@ptr
			} else {
				filter <- filter@ptr
			}
			if (proxy || is.null(extent)) {
				extent <- double()
			} else {
				extent <- as.vector(ext(extent))
			}
			p@ptr$read(x, layer, query, extent, filter, proxy)
			if (isTRUE(crs != "")) {
				crs(p) <- crs
			}
			if (proxy) {
				messages(p, "vect")
				pp <- methods::new("SpatVectorProxy")
				pp@ptr <- SpatVectorProxy$new()
				pp@ptr$v <- p@ptr
				return(pp)
			}
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
		methods::as(x, "SpatVector")
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
	if ((y[1] == "x") & (y[2] == "y")) return(FALSE)
	if ((x[1] == "eas") & (x[2] == "nor")) return(FALSE)
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
	return(FALSE)
}

setMethod("vect", signature(x="matrix"), 
	function(x, type="points", atts=NULL, crs="") {
		type <- tolower(type)
		type <- match.arg(tolower(type), c("points", "lines", "polygons"))
		stopifnot(NCOL(x) > 1)
		
		crs <- character_crs(crs, "vect")
		p <- methods::new("SpatVector")
		p@ptr <- SpatVector$new()
		crs(p) <- crs 

		nr <- nrow(x)
		if (nr == 0) {
			return(p)
		}

		if (ncol(x) == 2) { 
			lonlat <- .checkXYnames(colnames(x))
			if (type == "points") {
				p@ptr$setPointsXY(as.double(x[,1]), as.double(x[,2]))
			} else {
				p@ptr$setGeometry(type, rep(1, nr), rep(1, nr), x[,1], x[,2], rep(FALSE, nr))
			}
			if (lonlat && isTRUE(crs=="")) crs <- "+proj=longlat" 
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



setReplaceMethod("[", c("SpatVector", "ANY", "ANY"),
	function(x, i, j, value) {
		v <- values(x)
		v[i,j] <- value
		if (nrow(v) != nrow(x)) {
			error("[<-", "this would create an invalid SpatVector")
		}
		values(x) <- v
		x
	}
)

setReplaceMethod("[", c("SpatVector", "ANY", "missing"),
	function(x, i, j, value) {
		v <- values(x)
		if (inherits(value, "SpatVector")) {
			value <- values(value)
		}
		v[i,] <- value
		if (nrow(v) != nrow(x)) {
			error("[<-", "this would create an invalid SpatVector")
		}
		values(x) <- v
		x
	}
)

setReplaceMethod("[", c("SpatVector", "missing", "ANY"),
	function(x, i, j, value) {
		v <- values(x)
		if (inherits(value, "SpatVector")) {
			value <- values(value)
		}
		v[,j] <- value
		if (nrow(v) != nrow(x)) {
			error("[<-", "this would create an invalid SpatVector")
		}
		values(x) <- v
		x
	}
)


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
	function(x, geom=c("lon", "lat"), crs=NA) {
		if (!all(geom %in% names(x))) {
			error("vect", "the variable name(s) in argument `geom` are not in `x`")
		}
		crs <- character_crs(crs, "vect")
		if (length(geom) == 2) {
			geom <- match(geom[1:2], names(x))
			cls <- sapply(x[, geom], class)
			if (cls[1] == "integer") {
				x[,geom[1]] = as.numeric(x[,geom[1]])
			}
			if (cls[2] == "integer") {
				x[,geom[2]] = as.numeric(x[,geom[2]])
			}	
			p <- methods::new("SpatVector")
			p@ptr <- SpatVector$new()
			x <- .makeSpatDF(x)
			
			p@ptr$setPointsDF(x, geom-1, crs)
			messages(p, "vect")
			return(p)
		} else if (length(geom) == 1) {
			v <- vect(unlist(x[,geom]), crs=crs)
		} else {
			error("vect", "the length of 'geom' must be 1 or 2")
		}
		values(v) <- x
		v
	}
)


setMethod("vect", signature(x="list"), 
	function(x) {
		x <- svc(x)
		x <- x@ptr$append()
		v <- methods::new("SpatVector")
		v@ptr <- x 
		messages(v, "vect")
	}
)


setMethod("query", signature(x="SpatVectorProxy"), 
	function(x, start=1, n=nrow(x), vars=NULL, where=NULL, extent=NULL, filter=NULL, proxy=FALSE) {
		f <- x@ptr$v$source
		layer <- x@ptr$v$layer
		qy <- ""
		if (is.null(vars)) {
			vars <- "*"
		} else {
			vars <- na.omit(unique(vars))
			nms <- names(x)
			if (!all(vars %in% nms)) {
				error("query", "not all vars are variable names")
			} else if (length(vars) < length(nms))  {
				vars <- paste(vars, collapse=", ")
			}
		}

		if (!is.null(where)) {
			qy <- paste("SELECT", vars, "FROM", layer, "WHERE", where[1]) 
		}

		nr <- nrow(x)
		start <- start-1
		if (start > 0) {
			if (qy == "") {
				qy <- paste("SELECT", vars, "FROM", layer)
			} 
			if (n >= (nr-start)) {
				qy <- paste(qy, "OFFSET", start)
			} else {
				n <- min(n, nr-start)
				qy <- paste(qy, layer, "LIMIT", n, "OFFSET", start)
			}
		} else if (n < nr) {
			if (qy == "") {
				qy <- paste("SELECT", vars, "FROM", layer)
			} 
			n <- min(n, nr)
			qy <- paste(qy, "LIMIT", n)		
		}
		
		vect(f, layer, query=qy, extent=extent, filter=filter, crs="", proxy[1])
	}
)


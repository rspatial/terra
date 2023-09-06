
#setMethod("row.names", signature(x="SpatVector"),
#	function(x) {
#		1:nrow(x)
#	}
#)

setMethod("emptyGeoms", signature(x="SpatVector"),
	function(x) {
		x@cpp$nullGeoms() + 1
	}
)


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
		p@cpp <- SpatVector$new()
		messages(p, "vect")
		return(p)
	}
)

setMethod("vect", signature(x="SpatExtent"),
	function(x, crs="") {
		as.polygons(x, crs=crs)
	}
)

setMethod("vect", signature(x="SpatVectorCollection"),
	function(x) {
		vect(as.list(x))
	}
)


setMethod("vect", signature(x="character"),
	function(x, layer="", query="", extent=NULL, filter=NULL, crs="", proxy=FALSE, what="") {

		what <- trimws(tolower(what))
		if (what != "") what <- match.arg(trimws(tolower(what)), c("geoms", "attributes"))
		
		s <- substr(x[1], 1, 5)
		if (s %in% c("POINT", "MULTI", "LINES", "POLYG", "EMPTY")) {
			p <- methods::new("SpatVector")
#		if (all(grepl("\\(", x) & grepl("\\)", x))) {
			p@cpp <- SpatVector$new(gsub("\n", "", x))
			messages(p, "vect")
			crs(p, warn=FALSE) <- crs
			return(p)
		} 
		
		x <- x[1]
		nx <- try(normalizePath(x, mustWork=TRUE), silent=TRUE)
		if (!inherits(nx, "try-error")) { # skip html
			x <- nx
			if (grepl("\\.rds$", tolower(x))) {
				v <- unwrap(readRDS(x))
				if (!inherits(v, "SpatVector")) {
					error("vect", "the rds file does not store a SpatVector")
				}
				return(v)
			}
		} else if ((substr(x, 1, 4) == "http") & (grepl("\\.shp$", x) | grepl("\\.gpkg$", x))) {
			x <- paste0("/vsicurl/", x[1])
		}

		p <- methods::new("SpatVector")
		p@cpp <- SpatVector$new()
		proxy <- isTRUE(proxy)
			if ((what=="attributes") && proxy) {
			error("vect", "you cannot use 'what==attributes' when proxy=TRUE")
		}
		#if (proxy) query <- ""
		if (is.null(filter)) {
			filter <- SpatVector$new()
		} else {
			if (!inherits(filter, "SpatVector")) {
				error("vect", "'filter' must be a SpatVector")			
			}
			if (proxy) {
				error("vect", "you cannot use 'filter' when proxy=TRUE")
			}
			filter <- filter@cpp
		}
		if (is.null(extent)) {
			extent <- double()
		} else {
			extent <- as.vector(ext(extent))
		}
		p@cpp$read(x, layer, query, extent, filter, proxy, what)
		if (isTRUE(crs != "")) {
			crs(p, warn=FALSE) <- crs
		}
		if (proxy) {
			messages(p, "vect")
			pp <- methods::new("SpatVectorProxy")
			pp@cpp <- SpatVectorProxy$new()
			pp@cpp$v <- p@cpp
			return(pp)
		}
		p <- messages(p, "vect")
		if (what == "attributes") {
			p <- values(p)
		}
		p
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
	if ((y[1] == "x") && (y[2] == "y")) return(FALSE)
	if ((x[1] == "eas") && (x[2] == "nor")) return(FALSE)
	if ((x[1] == "lon") && (x[2] == "lat")) return(TRUE)
	if (grepl("lon", z[1]) && grepl("lat", z[2])) return(TRUE)

	if ((x[1] == "lat") && (x[2] == "lon")) {
		stop("vect", "longitude/latitude in the wrong order")
	} else if ((y[1] == "y") && (y[2] == "x")) {
		stop("vect", "x/y in the wrong order")
	} else if ((x[1] == "nor") && (x[2] == "eas")) {
		stop("vect", "easting/northing in the wrong order")
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
		p@cpp <- SpatVector$new()
		crs(p, warn=FALSE) <- crs

		nr <- nrow(x)
		if (nr == 0) {
			return(p)
		}
		nc <- ncol(x)
		if (nc == 2) {
			lonlat <- .checkXYnames(colnames(x))
			if (type == "points") {
				p@cpp$setPointsXY(as.double(x[,1]), as.double(x[,2]))
			} else {
				p@cpp$setGeometry(type, rep(1, nr), rep(1, nr), x[,1], x[,2], rep(FALSE, nr))
			}
			if (lonlat && isTRUE(crs=="")) crs <- "+proj=longlat"
		} else if (nc == 3) {
			p@cpp$setGeometry(type, x[,1], rep(1, nr), x[,2], x[,3], rep(FALSE, nr))
		} else if (nc == 4) {
			p@cpp$setGeometry(type, x[,1], x[,2], x[,3], x[,4], rep(FALSE, nr))
		} else if (nc == 5) {
			p@cpp$setGeometry(type, x[,1], x[,2], x[,3], x[,4], x[,5])
		} else {
			error("vect", "not an appropriate matrix (too many columns)")
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
	if (!(name %in% names(x))) {
		return(NULL)
		#error("$", paste(name, "is not a variable name in x"))
	}
	s <- .subset_cols(x, name, drop=TRUE)
	s[,1,drop=TRUE]
})


setMethod("[[", c("SpatVector", "numeric", "missing"),
function(x, i, j,drop=FALSE) {
	s <- .subset_cols(x, i, drop=TRUE)
	s[,,drop=drop]
})


setMethod("[[", c("SpatVector", "character", "missing"),
function(x, i, j, drop=FALSE) {
	if (!(any(i %in% names(x)))) {
		return(NULL)
	}
	s <- .subset_cols(x, i, drop=TRUE)
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


setReplaceMethod("[[", c("SpatVector", "character"),
	function(x, i, value) {

		x@cpp <- x@cpp$deepcopy()
		if (is.null(value)) {
			for (name in i) {
				if (name %in% names(x)) {
					x@cpp$remove_column(name)
				}
			}
			return(x);
		}

		if (inherits(value, "data.frame")) {
			if (ncol(value)	> 1) {
				warn("`[[<-`", "only using the first column")
			}
			value <- value[,1]
		} else if (inherits(value, "list")) {
			value <- unlist(value)
		}

		if (NCOL(value)	> 1) {
			warn("[[<-,SpatVector", "only using the first column")
			value <- value[,1]
		}
		name <- i[1]
		value <- rep(value, length.out=nrow(x))

		if (name %in% names(x)) {
			d <- values(x)
			if (all(is.na(value))) {
				#[] to keep type if NA is used
				d[[name]][] <- value
			} else {
				d[[name]] <- value			
			}
			values(x) <- d
		} else {
			if (inherits(value, "factor")) {
				v <- .makeSpatFactor(value)
				ok <- x@cpp$add_column_factor(v, name)
			} else if (inherits(value, "character")) {
				ok <- x@cpp$add_column_string(enc2utf8(value), name)
			} else if (inherits(value, "integer")) {
				# min long (should query what it is on the system?)
				value[is.na(value)] <- -2147483648
				ok <- x@cpp$add_column_long(value, name)
			} else if (inherits(value, "logical")) {
				v <- as.integer(value)
				v[is.na(v)] <- 2
				ok <- x@cpp$add_column_bool(v, name)
			} else if (inherits(value, "numeric")) {
				ok <- x@cpp$add_column_double(value, name)
			} else if (inherits(value, "Date")) {
				ok <- x@cpp$add_column_time(as.numeric(as.POSIXlt(value)), name, "days", "")
			} else if (inherits(value, "POSIXt")) {
				tz <- if (length(value) > 0) { attr(value[1], "tzone") } else { "" }
				if (is.null(tz)) tz <- ""
				ok <- x@cpp$add_column_time(as.numeric(value), name, "seconds", tz)
			} else {
				v <- try(as.character(value))
				if (!inherits(v, "try-error")) {
					ok <- x@cpp$add_column_string(enc2utf8(v), name)
				} else {
					ok <- FALSE
				}
			}
			if (!ok) {
				error("[[<-,SpatVector", "cannot add these values")
			}
		}
		x
	}
)


setReplaceMethod("[[", c("SpatVector", "numeric"),
	function(x, i, value) {
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
	function(x, geom=c("lon", "lat"), crs="", keepgeom=FALSE) {
		if (!all(geom %in% names(x))) {
			error("vect", "the variable name(s) in argument `geom` are not in `x`")
		}
		crs <- character_crs(crs, "vect")
		if (length(geom) == 2) {
			geom <- match(geom[1:2], names(x))
			if (inherits(x[,geom[1]], "integer")) {
				x[,geom[1]] = as.numeric(x[,geom[1]])
			}
			if (inherits(x[,geom[2]], "integer")) {
				x[,geom[2]] = as.numeric(x[,geom[2]])
			}
			p <- methods::new("SpatVector")
			p@cpp <- SpatVector$new()
			x <- .makeSpatDF(x)

			p@cpp$setPointsDF(x, geom-1, crs, keepgeom)
			messages(p, "vect")
			return(p)
		} else if (length(geom) == 1) {
			v <- vect(unlist(x[,geom]), crs=crs)
			if (!keepgeom) {
				x[[geom]] <- NULL
			}
		} else {
			error("vect", "the length of 'geom' must be 1 or 2")
		}
		values(v) <- x
		v
	}
)

setMethod("vect", signature(x="list"),
	function(x, type="points", crs="") {
		x <- lapply(x, function(i) {
			if (inherits(i, "SpatVector")) return(i)
			vect(i, type=type)
		})
		x <- svc(x)
		v <- methods::new("SpatVector")
		v@cpp <- x@cpp$append()
		if (crs != "") {
			crs(v) <- crs
		}
		messages(v, "vect")
	}
)




setMethod("query", signature(x="SpatVectorProxy"),
	function(x, start=1, n=nrow(x), vars=NULL, where=NULL, extent=NULL, filter=NULL, sql=NULL, what="") {
		f <- x@cpp$v$source
		slayer <- x@cpp$v$layer
		#1058
		layer <- paste0("\"", slayer, "\"")

		e <- x@cpp$v$read_extent
		if (is.null(extent)) {
			if (length(e) == 4) {
				extent = ext(e);
			}
		} else {
			if (length(e) == 4) {
				extent = intersect(ext(e), extent);
				if (is.null(extent)) {
					error("query", "extent does not intersect with x")
				}
			}
		}
		if (is.null(vars)) {
			vars <- "*"
		} else {
			vars <- stats::na.omit(unique(vars))
			nms <- names(x)
			if (!all(vars %in% nms)) {
				error("query", "not all vars are variable names")
			} else if (length(vars) < length(nms))  {
				vars <- paste(vars, collapse=", ")
			}
		}

		if (!is.null(sql)) {
			qy <- as.character(sql)
		
		} else {
			qy <- ""
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
		}
		
		if (qy != "") {
			if (x@cpp$v$read_query != "") {
				error("query", "A query was used to create 'x'; you can only subset it with extent or filter")
			}
		} else {
			layer <- slayer
		}

		p <- vect(f, layer, query=qy, extent=extent, filter=filter, crs="", FALSE, what=what)
		if (what == "attributes") {
			p <- values(p)
		}
		p
	}
)


vector_layers <- function(filename, delete="", return_error=FALSE) {
	p <- SpatVector$new()
	if (any(delete != "")) {
		delete <- trimws(delete)
		ok <- p$delete_layers(filename, delete, return_error[1])
		messages(p, "vector_layers")
		invisible(ok)
	} else {
		out <- p$layer_names(filename)
		messages(p, "vector_layers")
		out
	}
}


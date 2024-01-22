

setMethod("rasterizeGeom", signature(x="SpatVector", y="SpatRaster"),
	function(x, y, fun="count", unit="m", filename="", ...) {
		opt <- spatOptions(filename, ...)
		y@cpp <- y@cpp$rasterizeGeom(x@cpp, unit, fun, opt)
		messages(y, "rasterizeGeom")
	}
)

# now can use
# r@cpp = r@cpp$rasterizePoints(v@cpp, "mean", 1:nrow(v), NA, opt)

.set_names <- function(wopt, cnms, fun, nc) {
	if (is.null(wopt$names)) {
		if (is.null(cnms)) {
			if (nc == 1) {
				cnms <- fun
			} else {
				cnms <- paste0(fun, "_", 1:nc)
			}
		} else {
			cnms <- paste0(cnms, "_", fun)
		}
		wopt$names <- cnms
	}
	wopt
}

rasterize_points <- function(x, y, field, values, fun="last", background=NA, update=FALSE, filename="", overwrite=FALSE, wopt=list(), ...) {

	if (missing(fun)) fun <- "last"
	if (update && (!hasValues(y))) update <- FALSE
	nrx <- nrow(x)


	if (!is.data.frame(values)) {
		values <- as.data.frame(values)
	}
	if (nrow(values) == 1) {
		values <- sapply(values, function(x) rep_len(x, nrx))
		if (!is.data.frame(values)) { # dropped if nrx==1
			values <- as.data.frame(values)
		}
	} else if (nrow(values) != nrx) {
		error("rasterize", paste0("the number or rows in values is ", nrow(values), "\nThat does not match the number of points: ", nrx))
	}
	cnms <- colnames(values)
	nl <- ncol(values)
	r <- rast(y, nlyrs=nl)	
	levs <- list()
	has_levels <- FALSE
	for (i in 1:nl) {
		if (is.character(values[,i])) {
			f <- as.factor(values[,i])
			levs[[i]] <- levels(f)
			values[,i] <- as.integer(f) - 1
			has_levels <- TRUE
		} else if (is.factor(values[,i])) {
			f <- values[,i]
			levs[[i]] <- levels(f)
			values[,i] <- as.integer(f) - 1
			has_levels <- TRUE
		}
	}

	if (NCOL(values) == 1 && (!has_levels)) {
		txtfun <- .makeTextFun(fun)
		if (inherits(txtfun, "character")) {
			if (txtfun %in% c("first", "last", "pa", "sum", "mean", "count", "min", "max", "prod")) {	
				if (is.null(wopt$names)) {
					wopt$names <- txtfun
				}
				if (update) {
					ops <- spatOptions("", TRUE, wopt)	
				} else {
					ops <- spatOptions(filename, overwrite, wopt=wopt)				
				}
				narm <- isTRUE(list(...)$na.rm)
				r <- rast()
				r@cpp <- y@cpp$rasterizePointsXY(x[,1], x[,2], txtfun, values[[1]], narm, background, ops)
				messages(r)
				if (update) {
					r <- cover(r, y, filename=filename, overwrite=overwrite, wopt)
				}
				return(r)
			}
		}
	}
	if (inherits(fun, "character")) {
		if (fun == "first") {
			fun <- function(i, na.rm=FALSE) {
				if (na.rm) {
					i <- na.omit(i)
				}
				if (length(i) > 0) {
					i[i]
				} else {
					NA
				}
			}
		} else if (fun == "last") {
			fun <- function(i, na.rm=FALSE) {
				if (na.rm) {
					i <- na.omit(i)
				}
				if (length(i) > 0) {
					i[length(i)]
				} else {
					NA
				}
			}
		} else if (fun == "count") {
			fun <- function(i, na.rm=FALSE) {
				if (na.rm) {
					i <- na.omit(i)
				}
				length(i)
			}
		}
	}
	
	g <- cellFromXY(y, x)
	i <- which(!is.na(g))
	g <- g[i]
	if (length(g) == 0) {
		return(r)
	}
	values <- values[i, ,drop=FALSE]

	has_levels <- FALSE
	values <- aggregate(values, list(g), fun, ...)
	# allow for multiple fields
	#r[a[,1]] <- as.matrix(a[,-1])
	levs <- NULL
	if (is.null(wopt$names)) {
		fun <- .makeTextFun(fun)
		if (inherits(fun, "character")) {
			wopt <- .set_names(wopt, cnms, fun, NCOL(values))
		} else if (!is.null(cnms)) {
			wopt$names <- cnms
		}
	}

	values <- as.matrix(values)
	nl <- max(1, ncol(values)-1)
	r <- rast(r, nlyrs=nl)	

	if (!update) {
		if (has_levels) {
			levels(r) <- levs
		}
		b <- writeStart(r, filename=filename, sources=sources(y), overwrite=overwrite, wopt=wopt)
		filename  <- ""
	} else {
		b <- writeStart(r, "")
	}
	nc <- ncol(r)
	for (i in 1:b$n) {
		w <- matrix(background, nrow=b$nrows[i] * nc, ncol=nl)
		mincell <- cellFromRowCol(r, b$row[i], 1)
		maxcell <- cellFromRowCol(r, b$row[i] + b$nrows[i]-1, nc)
		vv <- values[values[,1] >= mincell & values[,1] <= maxcell, ,drop=FALSE]
		if (nrow(vv) > 0) {
			vv[,1] <- vv[,1] - (b$row[i] - 1) * nc
			w[vv[,1],] <- vv[,-1]
		}
		writeValues(r, w, b$row[i], b$nrows[i])
	}
	r <- writeStop(r)

	if (update) {
		r <- cover(r, y, filename=filename, overwrite=overwrite, wopt=wopt)
	}

	return (r)
}


setMethod("rasterize", signature(x="matrix", y="SpatRaster"),
	function(x, y, values=1, fun, ..., background=NA, update=FALSE, by=NULL, filename="", overwrite=FALSE, wopt=list()) {

		if (!is.null(by)) {
			by <- rep_len(by, nrow(x))
			values <- rep_len(values, nrow(x))

			x <- lapply(split(data.frame(x), by), as.matrix)
			values <- split(values, by)

			out <- rast(lapply(1:length(x), function(i) rasterize(x[[i]], y, values[[i]], fun, background=background, update=update)))
			names(out) <- unique(by)
			if (filename != "") {
				out <- writeRaster(out, filename, overwrite=overwrite, wopt=wopt)
			}
			return(out)
		}

		lonlat <- .checkXYnames(colnames(x))

		if (NCOL(values) <= 1) {
			values <- unlist(values)
			if (length(values) > nrow(x)) {
				error("rasterize", "length(values) > nrow(x)")
			}
			values=rep_len(values, nrow(x))
		} else {
			if (nrow(values) > nrow(x)) {
				error("rasterize", "nrow(values) > nrow(x)")
			}
			if (nrow(values) < nrow(x)) {
				i <- rep_len(1:nrow(values), nrow(x))
				values <- values[i, ]
			}
		}
		rasterize_points(x=x, y=y, field="", values=values, fun=fun, background=background, update=update, filename=filename, overwrite=overwrite, wopt=wopt, ...)
	}
)


setMethod("rasterize", signature(x="SpatVector", y="SpatRaster"),
	function(x, y, field="", fun, ..., background=NA, touches=FALSE, update=FALSE, cover=FALSE, by=NULL, filename="", overwrite=FALSE, wopt=list()) {

		if (!is.null(by)) {
			x <- split(x, by)
			uby <- sapply(x, function(i) i[[by]][1])			
			out <- rast(lapply(x, function(i) rasterize(i, y, field=field, fun, background=background, touches=touches, update=update, cover=cover, ...)))
			names(out) <- uby
			if (filename != "") {
				out <- writeRaster(out, filename, overwrite=overwrite, wopt=wopt)
			}
			return(out)
		}
		
		values <- 1
		if (!is.character(field)) {
			values <- as.numeric(field)
			field  <- ""
		} else if (is.na(field[1])) {
			values <- as.numeric(NA)
			field  <- ""
		} else if (is.null(field) || field[1] == "") {
			field <- ""
		} else if (!(field[1] %in% names(x))) {
			error("rasterize", paste(field, "is not a field in 'x'"))
		}

		g <- geomtype(x)
		if (grepl("points", g)) {
			nrx <- nrow(x)
			# also allow for multiple columns to multiple layers
			xy <- crds(x)
			if (field[1] == "") {
				values <- matrix(1, ncol=1, nrow=nrx)
			} else if (field[1] != "") {
				values <- x[, field, drop=TRUE]
				if (nrow(xy) != nrx) { # multi-points
					g <- geom(x)
					values <- values[g[,1], ,drop=FALSE]
				}
			}
			return(
				rasterize_points(x=xy, y=y, field=field, values=values, fun=fun, background=background, update=update, filename=filename, overwrite=overwrite, wopt=wopt, ...)
			)
		}

		opt <- spatOptions(filename, overwrite, wopt=wopt)
		pols <- grepl("polygons", g)

		if (cover[1] && pols) {
			y@cpp <- y@cpp$rasterize(x@cpp, "", 1, background, touches[1], "", TRUE, FALSE, TRUE, opt)
		} else {
			dots <- list(...)
			if (missing(fun)) {
				if (!is.null(dots$sum)) {
					# backward compatibility
					if (isTRUE(dots$sum)) fun <- "sum"
				} else {
					fun <- ""
				}
			}
			if (!inherits(fun, "character")) {
				fun <- .makeTextFun(fun)
				if (!inherits(fun, "character")) {
					error("rasterize", "'fun' must be 'min', 'max', 'mean', 'count', or 'sum'")
				}
			}
			if (fun != "") {
				fun <- tolower(fun)
				if (!(fun %in% c("sum", "mean", "min", "max", "count"))) {
					error("rasterize", "'fun' must be 'min', 'max' 'mean', 'count', or 'sum'")
				}
				if (fun == "count") {
					fun <- "sum"
					field <- ""
					values <- 1
				} else if (field != "") {
					if (fun == "min") {
						x <- sort(x[,field], field, TRUE)
						fun <- ""
					} else if (fun == "max") {
						x <- sort(x[,field], field, FALSE)
						fun <- ""
					}
				}
			}
			if ((field != "") && isTRUE(dots$na.rm)) {
				x <- x[!is.na(x[[field]]), ]
			}
			background <- as.numeric(background[1])
			if (fun == "sum") {
				xopt = spatOptions()
				y@cpp <- y@cpp$rasterize(x@cpp, field, values, background, touches[1], fun, FALSE, update[1], TRUE, xopt)
				messages(y, "rasterize")
				xopt = spatOptions()
				yy <- rast(y)
				yy@cpp <- y@cpp$rasterize(x@cpp, "", values, NA, touches[1], ""	, FALSE, update[1], TRUE, xopt)
				messages(yy, "rasterize")
				return(mask(y, yy, updatevalue=background, filename=filename, overwrite=overwrite, wopt=wopt))
			} else if (fun == "mean") {
				xopt = spatOptions()
				y@cpp <- y@cpp$rasterize(x@cpp, field, values, background, touches[1], "sum", FALSE, update[1], TRUE, xopt)
				messages(y, "rasterize")
				xopt = spatOptions()
				yy <- rast(y)
				yy@cpp <- y@cpp$rasterize(x@cpp, "", values, NA, touches[1], "sum", FALSE, update[1], TRUE, xopt)
				messages(yy, "rasterize")
				y <- y / yy
				if (filename != "") {
					y <- writeRaster(y, filename=filename, overwrite=overwrite, wopt=wopt)
				}
				return(y)
			} else {
				y@cpp <- y@cpp$rasterize(x@cpp, field, values, background, touches[1], fun, FALSE, update[1], TRUE, opt)
			}
		}	
		messages(y, "rasterize")
	}
)


setMethod("rasterize", signature(x="sf", y="SpatRaster"),
	function(x, y, ...) {
		x <- vect(x)
		rasterize(x, y, ...)
	}
)



rasterize_points <- function(x, y, field, values, fun="last", background=NA, update=FALSE, filename="", overwrite=FALSE, wopt=list(), ...) {

	if (update) {
		if (!hasValues(y)) {
			update <- FALSE
		} else {
			background <- NA 
		}
	} 

	nrx <- nrow(x)
	g <- geom(x, df=TRUE)
	# also allow for multiple columns to multiple layers
	if (field[1] != "") {
		values <- x[[field]]
	} else {
		if (!is.data.frame(values)) {
			values <- as.data.frame(values)
		}
		if (nrow(values) == 1) {
			values <- sapply(values, function(x) rep_len(x, nrx))
		} else {
			if (nrow(values) != nrx) {
				error("rasterize", "the number or rows in values does not match the number of points")
			}
		}
	}
	if (nrow(g) != nrx) { # multi-points
		values <- sapply(1:ncol(values), values[g[,1], i])
	}

	nl <- ncol(values)
	r <- rast(y, nlyrs=nl)
	values(r) <- background

	levs <- list()
	has_levels <- FALSE
	for (i in 1:nl) {
		if (is.character(values[,i])) {
			f <- as.factor(values[,i])
			levs[[i]] <- levels(f)
			values[,i] <- as.integer(f) - 1
			has_levels <- TRUE
		} 
	}
	
	g <- cellFromXY(y, as.matrix(g[, c("x", "y")]))
	i <- which(!is.na(g))
	g <- g[i]
	if (length(g) == 0) {
		return(r)
	}
	values <- values[i,]

	if (missing(fun)) fun <- "last"
	if (is.character(fun) && (fun %in% c("first", "last", "pa"))) {
		narm <- isTRUE(list(...)$na.rm)
		if (fun == "pa") {
			if (narm) {
				values <- aggregate(values, list(g), function(i) length(na.omit(i)))
				values[values < 1] <- background
			} else {
				values <- aggregate(values, list(g), function(i) 1)
			}
			has_levels <- FALSE
		} else if (fun == "first") {
			if (narm) {
				values <- aggregate(values, list(g), function(i) na.omit(i)[1])
			} else {
				values <- aggregate(values, list(g), function(i) i[1])
			}
		} else if (fun == "last") {
			if (narm) {
				values <- aggregate(values, list(g), function(i) rev(na.omit(i))[1])
			} else {
				values <- aggregate(values, list(g), function(i) rev(i)[1])
			}
		}
		r[values[,1]] <- as.matrix(values[,-1])

	} else {
		has_levels <- FALSE
		#a <- tapply(values, g, fun, ...)
		#b <- as.numeric(names(a))
		#r[b] <- as.vector(a)
		a <- aggregate(values, list(g), fun, ...)
		# allow for multiple fields
		r[a[,1]] <- as.matrix(a[,-1])
		levs <- NULL
	}

	if (update) {
		r <- cover(r, y)
	} else if (has_levels) {
		levels(r) <- levs
	}

	if (filename != "") {
		writeRaster(r, filename, overwrite=overwrite, wopt=wopt)
	}

	return (r)
}


setMethod("rasterize", signature(x="SpatVector", y="SpatRaster"), 
	function(x, y, field="", fun, ..., background=NA, touches=FALSE, update=FALSE, sum=FALSE, cover=FALSE, filename="", overwrite=FALSE, wopt=list()) {

		values <- 1
		if (is.null(field) || is.na(field) || (field == "")) {
			field <- ""
		} else if (!is.character(field)) {
			values <- as.numeric(field)
			field  <- ""
		} else {
			stopifnot(field %in% names(x))
		}

		g <- geomtype(x)
		if (grepl("points", g)) {
			r <- rasterize_points(x=x, y=y, field=field, values=values, fun=fun, background=background, update=update, filename=filename, overwrite=overwrite, wopt=wopt, ...) 
			return (r)
		}

		opt <- spatOptions(filename, overwrite, wopt=wopt)
		pols <- grepl("polygons", g)

		if (cover[1] && pols) {
			y@ptr <- y@ptr$rasterize(x@ptr, "", 1, background, touches[1], sum[1], TRUE, FALSE, TRUE, opt)
			y <- messages(y, "rasterize")
			return(y)
		}

		background <- as.numeric(background[1])
		y@ptr <- y@ptr$rasterize(x@ptr, field, values, background, touches[1], sum[1], FALSE, update[1], TRUE, opt)

		messages(y, "rasterize")
	}
)




# old_rasterize <- function(x, y, field, fun, background=NA, update=FALSE, touches=is.lines(x), cover=FALSE, filename="", ...) {

		# g <- geomtype(x)
		# if (grepl("points", g)) {
			# r <- rasterize_points(x=x, y=y, field=field, fun=fun, background=background, update=update, filename=filename, ...) 
			# return (r)
		# }

		# opt <- spatOptions(filename, ...)

		# if (cover[1] && grepl("polygons", g)) {
			# y@ptr <- y@ptr$rasterize1(x@ptr, "", 1, "", background, FALSE, touches[1], FALSE, TRUE, opt)
			# y <- messages(y, "rasterize")
			# return(y)
		# }

		# if (missing(field)) field <- 1:nrow(x)

		# domask <- FALSE
		# if (any(is.na(field))) {
			# if (length(field) == 1) {
				# y <- mask(y, x, inverse=FALSE, updatevalue=NA, touches=touches, filename=filename, ...)
				# return (y)
			# } else if (length(field) == nrow(x)) {
				# i <- is.na(field)
				# x <- x[!i, ]
				# if (update & !is.na(background)) {
					# xx <- x[i, ]
					# domask <- TRUE
					# fname <- filename
					# filename <- ""
				# }
			# } else {
				# error("rasterize", "NA detected in field, but length is not 1 or nrow(x)")
			# }
		# }


		# levs <- ""[0]
		# inverse <- FALSE # use "mask" for TRUE
		# background <- as.numeric(background[1])
		# #if (is.na(background)) background = 0/0 # NAN
		# if (is.numeric(field)) {
			# if ((length(field) == 1) && (ncol(x) > 1)) {
				# if (field > 0 && field <= ncol(x)) {
					# field <- names(x)[field]
				# } else {
					# error("rasterize", "field index outside the value range (1:ncol(x))")
				# } 
			# } else {
				# y@ptr <- y@ptr$rasterize1(x@ptr, "", field, levs, background, update[1], touches[1], inverse[1], FALSE, opt)
			# } 
		# }

		# if (is.character(field)) {
			# if (length(field) == 1) {
				# i <- which(field == names(x)) 
				# if (length(i) == 0) {
					# error("rasterize", paste0(field, " is not a valid fieldname"))
				# }
				# dtype <- datatype(x)[i]
				# if (dtype == "string") {
					# v <- x[[field]][,1]
					# f <- as.factor(v)
					# levs <- levels(f)
					# v <- as.integer(f) - 1
					# y@ptr <- y@ptr$rasterize1(x@ptr, field, v, levs, background, update[1], touches[1], inverse[1], FALSE, opt)
				# } else {
					
					# #y@ptr <- y@ptr$rasterize1(x@ptr, field, 0[0], levs, background, update[1], touches[1], inverse[1], opt)
					# ## for old gdal
					# y@ptr <- y@ptr$rasterize1(x@ptr, field, x[[field,drop=T]], levs, background, update[1], touches[1], inverse[1], FALSE, opt)
				# }
			# } else {
				# f     <- as.factor(field)
				# levs  <- levels(f)
				# field <- as.integer(f) - 1
				# y@ptr <- y@ptr$rasterize1(x@ptr, "value", field, levs, background, update[1], touches[1], inverse[1], FALSE, opt)
			# } 
		# } 

		# if (domask) {
			# y <- mask(y, xx, inverse=FALSE, updatevalue=NA, touches=touches, filename=fname, overwrite=overwrite, ...)
			# return (y)
		# }

		# messages(y, "rasterize")
	# }


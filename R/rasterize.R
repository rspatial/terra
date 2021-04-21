
rasterize_points <- function(x, y, field, values, fun="last", background, filename, update, ...) {

	
	if (update) {
		if (!hasValues(y)) {
			update <- FALSE
		} else {
			background <- NA 
		}
	} 
	r <- rast(y, nlyr=1)
	values(r) <- background

	g <- geom(x, df=TRUE)
	# also allow for multiple columns to multiple layers
	if (!missing(values)) {
		field <- rep_len(values, nrow(x))
	} else if (missing(field)) {
		field <- g[,1] # consider multi-point
	} else if (is.character(field)) {
		if ((length(field) == 1) && (field %in% names(x))) {
			names(r) <- field
			field <- as.vector(unlist(x[[field]]))
		} else if (length(field) != nrow(x)) {
			error("rasterize", "field length should be 1 or nrow(x)")
		}
	} else if (length(field) == 1) {
		if (field > 0 && field <= ncol(x)) {
			names(r) = names(x)[field]
			field <- as.vector(unlist(x[[field]]))
		} else {
			error("rasterize", "field index outside the value range (1:ncol(x))")
		}
	} else if (length(field) != nrow(x)) {
		error("rasterize", "length of field does not match the number of features")
	}

	levs <- NULL
	if (is.character(field)) {
		f <- as.factor(field)
		levs <- levels(f)
		field <- as.integer(f) - 1
	} 
	field <- field[g[,1]]

	g <- cellFromXY(y, as.matrix(g[, c("x", "y")]))
	i <- which(!is.na(g))
	g <- g[i]
	if (length(g) == 0) {
		return(r)
	}
	field <- field[i]

	if (missing(fun)) fun <- "last"
	if (is.character(fun)) {
		if (fun == "pa") {
			b <- unique(g)
			r[b] <- 1
		} else if (fun == "first") {
			r[rev(g)] <- rev(field)
		} else if (fun == "last") {
			r[g] <- field
		} else {
			error("rasterize", "unknown character function")
		}

	} else {
		a <- tapply(field, g, fun)
		b <- as.numeric(names(a))
		r[b] <- as.vector(a)
		levs <- NULL
	}

	if (update) {
		r <- cover(r, y)
	} else if (!is.null(levs)) {
		levels(r) <- levs
	}

	if (filename != "") {
		writeRaster(r, filename, ...)
	}

	return (r)
}


setMethod("rasterize", signature(x="SpatVector", y="SpatRaster"), 
	function(x, y, field="", values, fun, background=NA, touches=FALSE, update=FALSE, sum=FALSE, cover=FALSE, filename="", ...) {

		g <- geomtype(x)
		if (grepl("points", g)) {
			r <- rasterize_points(x=x, y=y, field=field, values=values, fun=fun, background=background, update=update, filename=filename, ...) 
			return (r)
		}

		opt <- spatOptions(filename, ...)
		pols <- grepl("polygons", g)
		#if (pols && touches) {
		#	first = rasterize(x, y, field=field, values=values, fun=fun, background=background, touches=FALSE, sum=sum, cover=cover, filename="", ...)
		#	ouf <- filename
		#	filename <- ""
		#}

		if (cover[1] && pols) {
			y@ptr <- y@ptr$rasterize(x@ptr, "", 1, background, touches[1], sum[1], TRUE, FALSE, opt)
			y <- messages(y, "rasterize")
			return(y)
		}

		if (is.null(field) || is.na(field)) field = ""
		stopifnot(length(field) == 1)
		stopifnot(is.numeric(values))
		if (field != "") {
			if (ncol(x) == 0) {
				error("rasterize", "there are no fields")
			}
			if (is.numeric(field)) {
				if (field > 0 && field <= ncol(x)) {
					field <- names(x)[field]
				} else {
					error("rasterize", "if `field` is numeric, it should be a single number in the range (1:ncol(x))")
				}
			} 
		} 
		background <- as.numeric(background[1])
		y@ptr <- y@ptr$rasterize(x@ptr, field, values, background, touches[1], sum[1], FALSE, update[1], opt)

		#if (pols && touches) {
		#	y <- cover(first, y, filename=filename, ...)
		#}
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


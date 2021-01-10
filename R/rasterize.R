
rasterize_points <- function(x=x, y=y, field=field, fun="last", background=background, update=update, filename=filename, overwrite=overwrite, wopt=wopt, ...) {

	if (update) {
		background <- NA 
	} 
	r <- rast(y, nlyr=1)
	values(r) <- background

	g <- geom(x, df=TRUE)
	# also allow for multiple columns to multiple layers
	if (missing(field)) {
		field <- g[,1] # consider multi-point
	} else if (is.character(field)) {
		if (length(field) == 1) {
			field <- as.vector(unlist(x[[field]]))
		} else if (length(field) != nrow(x)) {
			error("rasterize", "field length should be 1 or nrow(x)")
		}
	} else if (length(field) == 1) {
		if (field > 0 && field <= ncol(x)) {
			field <- as.vector(unlist(x[[field]]))
		} else {
			error("rasterize", "field index outside the value range (1:ncol(x))")
		}
	} else if (length(field) == nrow(x)) {
		field <- 1:nrow(x)
	} else  {
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
		a <- tapply(field, g, fun, ...)
		b <- as.numeric(names(a))
		r[b] <- as.vector(a)
	}

	if (update) {
		if (hasValues(y)) {
			r <- cover(r, y)
		}
	} else if (!is.null(levs)) {
		id <- (1:length(levs))-1
		cats(r, 1) <- data.frame(level=id, label=levs)
	}

	if (filename != "") {
		writeRaster(r, filename, overwrite=overwrite, wopt=wopt)
	}

	return (r)
}



setMethod("rasterize", signature(x="SpatVector", y="SpatRaster"), 
	function(x, y, field, fun, background=NA, update=FALSE, touches=is.lines(x), filename="", overwrite=FALSE, wopt=list(), ...) {

		if (grepl("points", geomtype(x))) {
			r <- rasterize_points(x=x, y=y, field=field, fun=fun, background=background, update=update, filename=filename, overwrite=overwrite, wopt=wopt, ...) 
			return (r)
		}

		if (missing(field)) field <- 1:nrow(x)

		domask <- FALSE
		if (any(is.na(field))) {
			if (length(field) == 1) {
				y <- mask(y, x, inverse=FALSE, updatevalue=NA, touches=touches, filename=filename, overwrite=overwrite, wopt=wopt, ...)
				return (y)
			} else if (length(field) == nrow(x)) {
				i <- is.na(field)
				x <- x[!i, ]
				if (update & !is.na(background)) {
					xx <- x[i, ]
					domask <- TRUE
					fname <- filename
					filename <- ""
				}
			} else {
				error("rasterize", "NA detected in field, but length is not 1 or nrow(x)")
			}
		}

		opt <- spatOptions(filename, overwrite, wopt)

		levs <- ""[0]
		inverse <- FALSE # use "mask" for TRUE
		background <- as.numeric(background[1])
		#if (is.na(background)) background = 0/0 # NAN
		if (is.character(field)) {
			if (length(field) == 1) {
				i <- which(field == names(x)) 
				if (length(i) == 0) {
					error("rasterize", paste0(field, " is not a valid fieldname"))
				}
				dtype <- datatype(x)[i]
				if (dtype == "string") {
					v <- x[[field]][,1]
					f <- as.factor(v)
					levs <- levels(f)
					v <- as.integer(f) - 1
					y@ptr <- y@ptr$rasterize(x@ptr, field, v, levs, background, update[1], touches[1], inverse[1], opt)
				} else {
					y@ptr <- y@ptr$rasterize(x@ptr, field, 0[0], levs, background, update[1], touches[1], inverse[1], opt)
				}
			} else {
				f     <- as.factor(field)
				levs  <- levels(f)
				field <- as.integer(f) - 1
				y@ptr <- y@ptr$rasterize(x@ptr, "value", field, levs, background, update[1], touches[1], inverse[1], opt)
			} 
		} else if (is.numeric(field)) {
			if (length(field) == 1 && (nrow(x) > 1)) {
				if (field > 0 && field <= ncol(x)) {
					v <- x[[field]][,1]
					if (inherits(v, "character")) {
						f <- as.factor(field)
						levs <- levels(f)
						v <- as.integer(f) - 1
					}
					y@ptr <- y@ptr$rasterize(x@ptr, field, v, levs, background, update[1], touches[1], inverse[1], opt)
				} else {
					error("rasterize", "field index outside the value range (1:ncol(x))")
				} 
			} else {
				y@ptr <- y@ptr$rasterize(x@ptr, "", field, levs, background, update[1], touches[1], inverse[1], opt)
			} 
		} else {
			error("rasterize", "field should be character or numeric")
		}

		if (domask) {
			y <- mask(y, xx, inverse=FALSE, updatevalue=NA, touches=touches, filename=fname, overwrite=overwrite, wopt=wopt, ...)
			return (y)
		}

		messages(y, "rasterize")
	}
)

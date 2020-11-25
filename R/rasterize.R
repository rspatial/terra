
rasterize_points <- function(x=x, y=y, field=field, fun="count", background=background, update=update, filename=filename, overwrite=overwrite, wopt=wopt, ...) {

	if (update) {
		background <- NA 
	} 
	r <- rast(y, nlyr=1)
	values(r) <- background

	g <- geom(x)
	
	# also allow for multiple columns to multiple layers
	if (missing(field)) {
		field <- g[,"id"] # consider multi-point
	} else if (is.character(field)) {
		field <- x[[field]]
		field <- field[g[,"id"], ,drop=FALSE]
	} else {
		if (length(field) == 1) {
			field <- rep(field, nrow(g))
		} else if (length(field) == nrow(x)) {
			field <- field[g[,"id"]]
		} else if (length(field) != nrow(g))  {
			stop("length of field does not match the number of features")
		}
	}
	g <- cellFromXY(y, as.matrix(g[, c("x", "y")]))
	
	if (missing(fun)) fun <- "pa"
	if (is.character(fun)) {
		if (fun == "pa") {
			b <- unique(g)
			r[b] <- 1
		} else if (fun == "first") {
			r[rev(g)] <- rev(field)
		} else if (fun == "last") {
			r[g] <- field
		} else {
			stop("unknown character function")
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
	}
	return (r)
}



setMethod("rasterize", signature(x="SpatVector", y="SpatRaster"), 
	function(x, y, field=1, fun, background=NA, update=FALSE, touches=is.lines(x), filename="", overwrite=FALSE, wopt=list(), ...) {

		if (geomtype(x) == "points") {
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
				stop("NA detected in field, but length is not 1 or nrow(x)")
			}
		}
		
		opt <- .runOptions(filename, overwrite, wopt)
					
		inverse=FALSE # use "mask" for TRUE
		background <- as.numeric(background[1])
		#if (is.na(background)) background = 0/0 # NAN
		if (is.character(field)) {
			y@ptr <- y@ptr$rasterize(x@ptr, field, 0, background, update[1], touches[1], inverse[1], opt)
		} else if (is.numeric(field)) {
			y@ptr <- y@ptr$rasterize(x@ptr, "", field, background, update[1], touches[1], inverse[1], opt)
		} else {
			stop("field should be character or numeric")
		}
		
		if (domask) {
			y <- mask(y, xx, inverse=FALSE, updatevalue=NA, touches=touches, filename=fname, overwrite=overwrite, wopt=wopt, ...)
			return (y)	
		}
		
		show_messages(y, "rasterize")
	}
)

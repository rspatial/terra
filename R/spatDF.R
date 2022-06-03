

.makeSpatDF <- function(d) {
	x <- methods::new("Rcpp_SpatDataFrame")
	nms <- colnames(d)
	for (i in 1:ncol(d)) {
		if (inherits(d[[i]], "factor")) {
			x$add_column_string(as.character(d[[i]]), nms[i])
		} else if (inherits(d[[i]], "character")) {
			x$add_column_string(enc2utf8(d[[i]]), nms[i])
		} else if (inherits(d[[i]], "integer")) {
			v <- d[[i]]
			# min long (should query what it is on the system?)
			v[is.na(v)] <- -2147483648
			x$add_column_long(v, nms[i])
		} else if (inherits(d[[i]], "logical")) {
			v <- as.integer(d[[i]])
			v[is.na(v)] <- 2
			x$add_column_bool(v, nms[i])
		} else if (inherits(d[[i]], "numeric")) {
			x$add_column_double(d[[i]], nms[i])
		} else if (inherits(d[[i]], "Date")) {
			x$add_column_time(as.numeric(as.POSIXlt(d[[i]])), nms[i], "days", "")
		} else if (inherits(d[[i]], "POSIXt")) {
			tz <- if (nrow(d) > 0) { attr(d[[i]][1], "tzone") } else { "" }
			x$add_column_time(as.numeric(d[[i]]), nms[i], "seconds", tz)
		} else {
			v <- try(as.character(d[[i]]))
			if (!inherits(v, "try-error")) {
				x$add_column_string(enc2utf8(v), nms[i])
			}
		}
	}
	x
}


.getSpatDF <- function(x, check.names = FALSE, stringsAsFactors=FALSE, ...) {
	d <- data.frame(x$values(), check.names=check.names, stringsAsFactors=stringsAsFactors, ...)
	if (ncol(d) == 0) return(d)
	
	s <- which(sapply(d, function(i) inherits(i, "character")))
	for (i in s) {
		d[[i]][d[[i]]=="NA"] <- NA
		Encoding(d[[i]]) <- "UTF-8"
	}
	ints <- which(x$itype == 1)
	for (i in ints) d[[i]] <- as.integer(d[[i]])
	bools <- which(x$itype == 3)
	for (i in bools) d[[i]] <- as.logical(d[[i]])
	times <- x$itype == 4
	if (any(times)) {
		steps <- x$get_timesteps()
		for (i in which(times)) {
			d[[i]] <- strptime("1970-01-01", "%Y-%m-%d", tz = "UTC") + d[[i]]
			if (steps[i] == "days") {
				d[[i]] <- as.Date(d[[i]])
			}
		}
	}
	d
}


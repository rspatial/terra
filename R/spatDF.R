

.makeSpatDF <- function(d) {
	x <- methods::new("Rcpp_SpatDataFrame")
	cls <- sapply(d, class)
	nms <- colnames(d)
	for (i in 1:length(cls)) { 
		if (cls[i] == "factor") {
			x$add_column_string(as.character(d[[i]]), nms[i])
		} else if (cls[i] == "character") {
			x$add_column_string(d[[i]], nms[i])
		} else if (cls[i] == "integer") {
			v <- d[[i]]
			# min long
			v[is.na(v)] <- -2147483648
			x$add_column_long(v, nms[i])
		} else {
			v <- as.numeric(d[[i]])
			x$add_column_double(v, nms[i])
		}
	}
	x
}


.getSpatDF <- function(x, check.names = FALSE, stringsAsFactors=FALSE, ...) {
	d <- data.frame(x$values(), check.names=check.names, stringsAsFactors=stringsAsFactors, ...)
	d[d=="NA"] <- NA
	s <- which(sapply(d, class) == "character")
	for (i in s) Encoding(d[[i]]) <- "UTF-8"
	ints <- which(x$itype == 1)
	for (i in ints) d[[i]] <- as.integer(d[[i]])
	d
}




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
			x$add_column_long(d[[i]], nms[i])
		} else {
			v <- as.numeric(d[[i]])
			x$add_column_double(v, nms[i])
		}
	}
	x
}


.getSpatDF <- function(x) {
	d <- data.frame(x$values(), check.names = FALSE, stringsAsFactors=FALSE)
	d[d=="NA"] <- NA
	d
}



setMethod("meta", signature(x="SpatRaster"),
	function(x, layers=FALSE) {
		f <- function(i) {
			if (length(i) == 0) {
				matrix(ncol=2, nrow=0)
			} else {
				matrix(unlist(regmatches(i, regexpr("=", i), invert=TRUE)), ncol=2, byrow=TRUE)
			}
		}
		lapply(x@cpp$metadata(layers), f)
	}
)


setMethod("metags", signature(x="SpatRaster"),
	function(x, layer=NULL, name=NULL) {
		if (!is.null(layer)) {
			v <- x@cpp$getLyrTags(layer-1)
		} else {
			v <- x@cpp$getTags()
		}
		m <- matrix(v, ncol=2, byrow=TRUE, dimnames = list(NULL, c("name", "value")))
		out <- m[,2]
		names(out) <- m[,1]
		if (!is.null(name)) {
			out <- out[name]
		} 
		out
	}
)


setMethod("metags<-", signature(x="SpatRaster"),
	function(x, ..., layer=NULL, value) {
		if (is.null(value)) {
			if (!is.null(layer)) {
				value <- matrix(x@cpp$getLyrTags(layer-1), ncol=2, byrow=TRUE)
			} else {
				value <- matrix(x@cpp$getTags(), ncol=2, byrow=TRUE)
			}
			value[,2] <- ""
			value[is.na(value)] <- ""
		} else if (NCOL(value) == 1) {
			if (!is.null(names(value)) && (!any(grepl("=", value)))) {
				value <- cbind(names(value), value)	
			} else {	
				value <- strsplit(value, "=")
				i <- sapply(value, length) == 1
				if (length(i) > 0) {
					j <- which(i)
					for (i in j) value[[i]] <- c(value[[i]], "")
				}
				i <- sapply(value, length) == 2
				value <- do.call(rbind, value[i])
			}
		} else if (NCOL(value) != 2) {
			error("metags<-", "expecting a vector with 'name=value' or a two column matrix")
		}
		value[is.na(value[,2]), 2] <- ""
		value <- na.omit(value)
		x <- deepcopy(x)
		if (NROW(value) > 0) {
			if (!is.null(layer)) {
				sapply(1:nrow(value), function(i) {
						x@cpp$addLyrTags(layer-1, value[i,1], value[i,2])
					})
			} else {
				sapply(1:nrow(value), function(i) {
						x@cpp$addTag(value[i,1], value[i,2])
					})
			}
		}
		x
	}
)


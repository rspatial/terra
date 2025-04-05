
setMethod("meta", signature(x="SpatRaster"),
	function(x, layers=FALSE) {
		f <- function(i) {
			if (length(i) == 0) {
				matrix(ncol=2, nrow=0)
			} else {
				matrix(unlist(regmatches(i, regexpr("=", i), invert=TRUE)), ncol=2, byrow=TRUE)
			}
		}
		lapply(x@pntr$metadata(layers), f)
	}
)


setMethod("metags", signature(x="SpatRaster"),
	function(x, layer=NULL, name=NULL) {
		if (!is.null(layer)) {
			if (is.character(layer)) layer = match(layer, names(x))
			v <- x@pntr$getLyrTags(layer-1)
			out <- matrix(v, ncol=3, byrow=TRUE, dimnames = list(NULL, c("layer", "name", "value")))
			out <- data.frame(out)
			out$layer <- as.numeric(out$layer) + 1
			if (!is.null(name)) {
				out <- out[out$name == name, , drop=FALSE]
			} 
		} else {
			out <- do.call(cbind, x@pntr$getTags())
			if (is.null(out)) return(out)
			colnames(out) <- c("domain", "name", "value")
			out <- data.frame(out)
			if (!is.null(name)) {
				out <- out[out$name == name, ]
			} 
			out <- out[, c(2,3,1)]
		}
		out
	}
)


parse_tags <- function(value, domain) {
	if (NCOL(value) == 1) {
		if (!is.null(names(value)) && (!any(grepl("=", value)))) {
			value <- cbind(names(value), value)	
		} else {	
			
			val <- strsplit(value, "=")
			i <- sapply(val, length) == 1
			if (sum(i) > 0) {
				j <- which(i)
				for (i in j) val[[i]] <- c(val[[i]], "")
			}
			i <- sapply(val, length) == 2
			val <- do.call(rbind, val[i])

			dom <- strsplit(val[,1], ":")
			n <- sapply(dom, length) 
			i <- n == 1
			if (sum(i) > 0) {
				for (i in which(i)) dom[[i]] <- c("", dom[[i]])
			}
			n <- sapply(dom, length) 
			i <- n > 2
			if (sum(i) > 0) {
				for (i in which(i)) dom[[i]] <- c(dom[[i]][1], paste(dom[[i]][2:length(dom[[i]])], collapse=":"))
			}
			i <- sapply(dom, length) == 2
			dom <- do.call(rbind, dom[i])
			value <- cbind(dom, val[,2])
		}
	} else if (NCOL(value) == 3) {
		value <- value[, c(3,1,2), drop=FALSE]
	} else if (NCOL(value) > 3) {
		error("metags<-", "expecting a vector with 'name=value' or a two/three column matrix")
	} 
	if (NCOL(value) == 2) value <- cbind(domain, value) 
	value[is.na(value[,2]), 2] <- ""
	na.omit(value)
}



setMethod("metags<-", signature(x="SpatRaster"),
	function(x, ..., layer=NULL, domain="", value) {
		if (is.null(value)) {
			if (!is.null(layer)) {
				if (is.character(layer)) layer = match(layer, names(x))
				value <- metags(x, layer)
			} else {
				value <- metags(x)
			}
			value[,2] <- ""
			#value[is.na(value)] <- ""
		} else {
			value <- parse_tags(value, domain)
		}
		x <- deepcopy(x)
		if (NROW(value) > 0) {
			if (!is.null(layer)) {
				if (is.character(layer)) layer = match(layer, names(x))
				x@pntr$addLyrTags(layer-1, value[,2], value[,3])
			} else {
				sapply(1:nrow(value), function(i) {
						x@pntr$addTag(value[i,2], value[i,3], value[i,1])
					})
			}
		}
		x
	}
)



setMethod("metags", signature(x="SpatRasterDataset"),
	function(x, dataset=NULL, name=NULL) {
		if (!is.null(dataset)) {
			if (is.character(dataset)) layer = match(dataset, names(x))		
			return(metags(x[[dataset]], name=name))
		} else {
			v <- x@pntr$getTags()
			m <- matrix(v, ncol=2, byrow=TRUE, dimnames = list(NULL, c("name", "value")))
			out <- m[,2]
			names(out) <- m[,1]
			if (!is.null(name)) {
				out <- out[name]
			} 
		}
		out
	}
)


setMethod("metags<-", signature(x="SpatRasterDataset"),
	function(x, ..., dataset=NULL, value) {
		if (is.null(value)) {
			if (!is.null(dataset)) {
				if (is.character(dataset)) layer = match(dataset, names(x))		
				value <- matrix(x[[dataset]]@pntr$getTags(), ncol=2, byrow=TRUE)
			} else {
				value <- matrix(x@pntr$getTags(), ncol=2, byrow=TRUE)
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
		x@pntr <- x@pntr$deepcopy()
		if (NROW(value) > 0) {
			if (!is.null(dataset)) {
				if (is.character(dataset)) layer = match(dataset, names(x))		
				x[[dataset]]@pntr$addTag(value[,1], value[,2])
			} else {
				sapply(1:nrow(value), function(i) {
						x@pntr$addTag(value[i,1], value[i,2])
					})
			}
		}
		x
	}
)


setMethod("metags", signature(x="SpatRasterCollection"),
	function(x, dataset=NULL, name=NULL) {
		if (!is.null(dataset)) {
			return(metags(x[[dataset]], name=name))
		} else {
			v <- x@pntr$getTags()
			m <- matrix(v, ncol=2, byrow=TRUE, dimnames = list(NULL, c("name", "value")))
			out <- m[,2]
			names(out) <- m[,1]
			if (!is.null(name)) {
				out <- out[name]
			} 
		}
		out
	}
)


setMethod("metags<-", signature(x="SpatRasterCollection"),
	function(x, ..., dataset=NULL, value) {
		if (is.null(value)) {
			if (!is.null(dataset)) {
				value <- matrix(x[[dataset]]@pntr$getTags(), ncol=2, byrow=TRUE)
			} else {
				value <- matrix(x@pntr$getTags(), ncol=2, byrow=TRUE)
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
		x@pntr <- x@pntr$deepcopy()
		if (NROW(value) > 0) {
			if (!is.null(dataset)) {
				x[[dataset]]@pntr$addTag(value[,1], value[,2])
			} else {
				sapply(1:nrow(value), function(i) {
						x@pntr$addTag(value[i,1], value[i,2])
					})
			}
		}
		x
	}
)


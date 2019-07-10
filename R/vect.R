

if (!isGeneric("vect") ) { setGeneric("vect", function(x, ...) standardGeneric("vect")) }

setMethod("vect", signature(x="missing"), 
	function(...) {
		p <- methods::new("SpatVector")
		p@ptr <- SpatVector$new()
		show_messages(p)
		return(p)
	}
)

setMethod("vect", signature(x="character"), 
	function(x, ...) {
		p <- methods::new("SpatVector")
		p@ptr <- SpatVector$new()
		x <- normalizePath(x)
		p@ptr$read(x)
		show_messages(p)
	}
)


setMethod("vect", signature(x="Spatial"), 
	function(x, ...) {
		methods::as(x, "SpatVector")
	}
)


setMethod("vect", signature(x="matrix"), 
	function(x, type="points", atts=NULL, crs=NA, ...) {
		type <- tolower(type)
		stopifnot(type %in% c("points", "lines", "polygons"))
		
		p <- methods::new("SpatVector")
		p@ptr <- SpatVector$new()
		nr <- nrow(x)
		if (nr == 0) {
			return(p)
		}
		
		if (ncol(x) == 2) {
			p@ptr$setGeometry(type, 1:nr, rep(1, nr), x[,1], x[,2], rep(FALSE, nr))
		} else if (ncol(x) == 4) {
			p@ptr$setGeometry(type, x[,1], x[,2], x[,3], x[,4], rep(FALSE, nr))		
		} else if (ncol(x) == 5) {
			p@ptr$setGeometry(type, x[,1], x[,2], x[,3], x[,4], x[,5])
		} else {
			stop("not an appropriate matrix")
		}
		if (!is.null(atts)) {
			types <- sapply(atts, class)
			nms <- colnames(atts)
			for (i in 1:ncol(atts)) {
				if (types[i] == "numeric") {
					p@ptr$add_column_double(atts[[i]], nms[i])
				} else if (types[i] == "integer") {
					p@ptr$add_column_long(atts[[i]], nms[i])
				} else if (types[i] == "character") {
					p@ptr$add_column_string(atts[[i]], nms[i])
				} else {
					att <- as.character(atts[[i]])
					p@ptr$add_column_string(att, nms[i])
				}
			}
		}
		if (!is.na(crs)) {
			p@ptr$crs <- ifelse(is.na(crs), "", as.character(crs))
		}
		show_messages(p)
	}
)


setMethod("vect", signature(x="data.frame"), 
	function(x, type="points", atts=NULL, crs=NA, ...) {
		x <- as.matrix(x)
		vect(x, type=type, atts=atts, crs=crs, ...)
	}
)

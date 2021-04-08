# Author: Robert J. Hijmans
# License GPL v3


setMethod("match", signature(x="SpatRaster"),
	function(x, table, nomatch=NA, incomparables=NULL) {
		table <- unique(table)
		if (is.character(table)) {
			if (!is.factor(x)) error("match", "table has character values")
			table <- stats::na.omit(match(table, levels(x)[[1]]))
			if (length(table) == 0) {
				return (x * 0)
			} 
		}
		app(x, function(i) match(i, table, nomatch, incomparables))
	}
)


setMethod("%in%", signature(x="SpatRaster"),
	function(x, table) {
		opt <- spatOptions("", FALSE, list())
		table <- unique(table)
		x@ptr <- x@ptr$is_in(table, opt)
		messages(x, "%in%")
	}
)

#setMethod("%in%", signature(x="SpatRaster", table="ANY"),
#	function(x, table) {
#		out <- rast(x)
#		readStart(x)
#		on.exit(readStop(x))
#		nc <- ncol(out)
#		b <- writeStart(out, filename="", overwrite=FALSE, wopt=list())
#		for (i in 1:b$n) {
#			v <- readValues(x, b$row[i], b$nrows[i], 1, nc, TRUE)
#			v <- v %in% table
#			writeValues(out, v, b$row[i], b$nrows[i])
#		}
#		writeStop(out)
#	}
#)



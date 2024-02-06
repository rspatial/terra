# Author: Robert J. Hijmans
# License GPL v3

getCatIDs <- function(x, table, sender="%in%") {
	if (!is.factor(x)) {
		error(sender, "Can only match character values if x is categorical")
	}
	if (nlyr(x) != 1) {
		error(sender, "matching with character values is only supported for single layer SpatRaster")
	}
	d <- cats(x)[[1]]
	m <- stats::na.omit(match(table, d[,2]))
	d[m,1]
}


setMethod("match", signature(x="SpatRaster"),
	function(x, table, nomatch=NA, incomparables=NULL) {
		table <- unique(table)
		if (is.character(table)) {
			table <- getCatIDs(x, table, sender="match")
			if (length(table) == 0) {
				return(as.logical(x*0))
			}
		}
		app(x, function(i) match(i, table, nomatch, incomparables))
	}
)


setMethod("%in%", signature(x="SpatRaster"),
	function(x, table) {

		table <- unique(table)
		if (is.character(table)) {
			table <- getCatIDs(x, table)
			if (length(table) == 0) {
				return(as.logical(x*0))
			}
		}
		opt <- spatOptions("", FALSE, list())
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
#		b <- writeStart(out, filename="", overwrite=FALSE, sources=sources(x), wopt=list())
#		for (i in 1:b$n) {
#			v <- readValues(x, b$row[i], b$nrows[i], 1, nc, TRUE)
#			v <- v %in% table
#			writeValues(out, v, b$row[i], b$nrows[i])
#		}
#		writeStop(out)
#	}
#)



# Author: Robert J. Hijmans
# License GPL v3

getFactTable <- function(x, table, sender="%in%") {
	if (!is.factor(x)) {
		error(sender, "Can only match character values of x is categorical")
	}
	if (nlyr(x) != 1) {
		error(sender, "matching with character values is only supported for single layer SpatRaster")
	}
	d <- cats(x)[[1]]
	levs <- levels(x)[[1]]
	m <- na.omit(match(table, levs))
	if (length(m) == 0) {
		return(as.logical(x*0))
	}
	d[m,1]
}

setMethod("match", signature(x="SpatRaster"),
	function(x, table, nomatch=NA, incomparables=NULL) {
		table <- unique(table)
		if (is.character(table)) {
			table <- getFactTable(x, table, sender="match")
		}
		app(x, function(i) match(i, table, nomatch, incomparables))
	}
)


setMethod("%in%", signature(x="SpatRaster"),
	function(x, table) {

		table <- unique(table)
		if (is.character(table)) {
			table <- getFactTable(x, table)
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
#		b <- writeStart(out, filename="", overwrite=FALSE, wopt=list())
#		for (i in 1:b$n) {
#			v <- readValues(x, b$row[i], b$nrows[i], 1, nc, TRUE)
#			v <- v %in% table
#			writeValues(out, v, b$row[i], b$nrows[i])
#		}
#		writeStop(out)
#	}
#)



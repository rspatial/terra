# Author: Robert J. Hijmans
# License GPL v3


setMethod("match", signature(x="SpatRaster", table="ANY", nomatch="ANY", incomparables="ANY"),
	function(x, table, nomatch, incomparables) {
#		app(x, function(i) match(i, table, nomatch, incomparables))
		capp(x, function(i) match(i, table, nomatch, incomparables))
	}
)


setMethod("%in%", signature(x="SpatRaster", table="ANY"),
	function(x, table) {
		opt <- .runOptions("", FALSE, list())
		table <- unique(table)
		x@ptr <- x@ptr$is_in(table, opt)
		show_messages(x, "%in%")
	}
)

#setMethod("%in%", signature(x="SpatRaster", table="ANY"),
#	function(x, table) {
#		out <- rast(x)
#		if (!readStart(x)) { stop(x@ptr$messages$getError()) }
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



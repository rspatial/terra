# Author: Robert J. Hijmans
# Date : January 2009 - December 2011
# Version 1.0
# Licence GPL v3




#.getClicks <- function(...) {
#	res <- list()
#	while(TRUE) {
#		loc <- graphics::locator(1, ...)
#		if (is.null(loc)) break
#		res <- c(res, loc)
#	}
#	matrix(res, ncol=2, byrow=TRUE)
#}



.getCellFromClick <- function(x, n, type, id, ...) {
	#loc <- graphics::locator(n, type, ...)
	#xyCoords <- cbind(x=loc$x, y=loc$y)
	xyCoords <- RS_locator(n, type, ...)
	if (id) {
		text(xyCoords, labels=1:n)
	}
	cells <- cellFromXY(x, xyCoords)
	cells <- unique(stats::na.omit(cells))
	if (length(cells) == 0 ) { 
		error("click", "no valid cells selected") 
	}
	cells
}



setMethod("click", signature(x="missing"), 
	function(x, n=1, id=FALSE, type="p", ...) {
		#loc <- graphics::locator(n, type, ...)
		#cbind(x=loc$x, y=loc$y)
		RS_locator(n, type=type, id=id, ...)
	}
)


setMethod("click", signature(x="SpatRaster"), 
	function(x, n=1000, id=FALSE, xy=FALSE, cell=FALSE, type="p", show=TRUE, ...) {
	n <- max(n, 1)
	xyCoords <- RS_locator(n, type=type, id=id, ...)
	if (NROW(xyCoords) == 0) return(invisible(NULL))
	values <- NULL
	cells <- stats::na.omit(cellFromXY(x, xyCoords))
	if (length(cells) == 0) break
	value <- x[cells]
	if (cell) {
		value <- data.frame(cell=cells, value)
	}
	if (xy) { 
		xyCoords <- xyFromCell(x, cells)
		colnames(xyCoords) <- c("x", "y")
		value <- data.frame(xyCoords, value)
	} 
	if (show) {
		print(value)
		utils::flush.console()
	}
	if (is.null(dim(value))) { 
		value <- matrix(value)
		colnames(value) <- names(x)
	}
	values <- rbind(values, value)
	if (show) {
		invisible(values)
	} else {
		values
	}
})



setMethod("click", signature(x="SpatVector"), 
	function(x, n=1, id=FALSE, type="p", ...) {
		#loc <- graphics::locator(n, type, ...)
		#xy <- vect(cbind(x=loc$x, y=loc$y))
		xy <- vect(RS_locator(n, type=type, id=id, ...))
		e <- extract(xy, x)
		if (!id) {
			e[,-1]
		} else {
			colnames(e)[1] <- "ID"
			e
		}
	}
)



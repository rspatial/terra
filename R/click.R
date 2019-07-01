# Author: Robert J. Hijmans
# Date : January 2009 - December 2011
# Version 1.0
# Licence GPL v3




.getClicks <- function(...) {
	res <- list()
	while(TRUE) {
		loc <- graphics::locator(1, ...)
		if (is.null(loc)) break
		res <- c(res, loc)
	}
	matrix(res, ncol=2, byrow=TRUE)
}



.getCellFromClick <- function(x, n, type, id, ...) {
	loc <- graphics::locator(n, type, ...)
	xyCoords <- cbind(x=loc$x, y=loc$y)
	if (id) {
		text(xyCoords, labels=1:n)
	}
	cells <- cellFromXY(x, xyCoords)
	cells <- unique(stats::na.omit(cells))
	if (length(cells) == 0 ) { 
		stop('no valid cells selected') 
	}
	cells
}



setMethod('click', signature(x='missing'), 
	function(x, n=1, type="n", ...) {
		loc <- graphics::locator(n, type, ...)
		cbind(x=loc$x, y=loc$y)
	}
)

	
setMethod('click', signature(x='SpatRaster'), 
	function(x, n=Inf, id=FALSE, xy=FALSE, cell=FALSE, type="n", show=TRUE, ...) {
	values <- NULL
	i <- 0
	n <- max(n, 1)
	while (i < n) {
		i <- i + 1
		loc <- graphics::locator(1, type, ...)
		if (is.null(loc)) break
		xyCoords <- cbind(x=loc$x, y=loc$y)
		if (id) { 
			text(xyCoords, labels=i) 
		}
		cells <- stats::na.omit(cellFromXY(x, xyCoords))
		if (length(cells) == 0) break
		
		value <- x[cells]
		if (cell) {
			value <- data.frame(cell=cells, value)
		}
		if (xy) { 
			xyCoords <- xyFromCell(x, cells)
			colnames(xyCoords) <- c('x', 'y')
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
	}
	if (show) {
		invisible(values)
	} else {
		values
	}
})




setMethod('click', signature(x='SpatVector'), 
	function(x, ...) {
		stop("not implemented yet")
}		
)


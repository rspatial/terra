# Author: Robert J. Hijmans
# Date :  December 2020
# Version 1.0
# License GPL v3

.inset <- function(x, e, loc="", scale=0.2, background="white", border="black", box=NULL, pbx, ...) {

	usr <- par("usr")
	if (missing(e)) {
		e <- ext(usr)
		r <- diff(e[1:2]) / diff(e[3:4])
		e[2] <- e[1] + scale * diff(e[1:2])
		e[3] <- e[4] - scale * diff(e[3:4]) * r
	}

	offset <- 0.9
	#offset <- max(0.1, min(1, offset))
	scale  <- offset * min(e / ext(x))
	
	y  <- rescale(x, scale)
	ey <- ext(y)
	xy <- c(mean(ey[1:2]), mean(ey[3:4]))
	xybox <- c(mean(e[1:2]), mean(e[3:4]))
	dx <- xybox[1] - xy[1]
	dy <- xybox[2] - xy[2]
	y  <- shift(y, dx, dy)
	if (!is.null(box)) {
		ex <- ext(x)
		box  <- rescale(as.polygons(box), scale, ex[1]+diff(ex[1:2])/2, ex[3]+diff(ex[3:4])/2)
		box <- shift(box, dx, dy)
	}
	
	if (loc != "") {
		stopifnot(loc %in% c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center"))
		ex <- ext(x)
		if (grepl("top", loc)) {
			dy <- usr[4] - e[4]
		} else if (grepl("bottom", loc)) {
			dy <- usr[3] - e[3]
		} else {
			dy <- (usr[3] + diff(usr[3:4])/2) - (e[3] + diff(e[3:4])/2)
		}
		if (grepl("left", loc)) {
			dx <- usr[1] - e[1]		
		} else if (grepl("right", loc)) {
			dx <- usr[2] - e[2]		
		} else {
			dx <- (usr[1] + diff(usr[1:2])/2) - (e[1] + diff(e[1:2])/2)
		}
		y <- shift(y, dx, dy)
		e <- shift(e, dx, dy)
		if (!is.null(box)) {
			box <- shift(box, dx, dy)
		}
	}
	if (!is.na(background)) {
		polys(as.polygons(e), col=background, lty=0)
	}

	plot(y, ..., axes=FALSE, legend=FALSE, add=TRUE)
	if (!is.na(border)) {
		lines(e, col=border)
	}
	if (!is.null(box)) {
		if (missing(pbx) || !is.list(pbx)) {
			pbx <- list()
		}
		pbx$x <- box
		do.call(lines, pbx)
	}
	invisible(y)
}

setMethod("inset", signature(x="SpatVector"), 
	function(x, e, loc="", scale=0.2, background="white", border="black", box=NULL, pbx, ...) {
		.inset(x, e, loc=loc, scale=scale, background=background, border=border, box=box, pbx=pbx, ...)
	}
)


setMethod("inset", signature(x="SpatRaster"), 
	function(x, e, loc="", scale=0.2, background="white", border="black", box=NULL, pbx, ...) {
		.inset(x, e, loc=loc, scale=scale, background=background, border=border, box=box, pbx=pbx, ...)
	}
)

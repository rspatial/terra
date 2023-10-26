# Author: Robert J. Hijmans
# Date :  December 2020
# Version 1.0
# License GPL v3


setMethod("inext", signature(x="SpatVector"),
	function(x, e, y=NULL, gap=0) {
	# the area used for scaling
		gap <- rep_len(gap, 2)
		e <- as.vector(e) + c(gap[1], -gap[1], gap[2], -gap[2])
		stopifnot((e[2] > e[1]) && (e[4] > e[3]))
		ex <- ext(x)
		x <- shift(x, e[1] - ex[1], e[3] - ex[3])
		ve <- ext(x)
		fx <- (e[2] - e[1]) / (ve[2] - ve[1])
		fy <- (e[4] - e[3]) / (ve[4] - ve[3])

		if (!is.null(y)) {
			y <- shift(y, e[1] - ex[1], e[3] - ex[3])
			rescale(y, fx=fx, fy=fy, e[1], e[3])
		} else {
			rescale(x, fx=fx, fy=fy, e[1], e[3])
		}
	}
)


.inset <- function(x, e, loc="", scale=0.2, background="white", perimeter="black", pper, box=NULL, pbox, xpd=NA, ...) {

	usr <- unlist(get.clip()[1:4])
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
		box <- rescale(as.polygons(box), scale, x0=ex[1]+diff(ex[1:2])/2, y0=ex[3]+diff(ex[3:4])/2)
		box <- shift(box, dx, dy)
	}

	if ((loc != "") && (loc != "topleft")) {
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
		polys(as.polygons(e), col=background, lty=0, xpd=xpd)
	}

	plot(y, ..., axes=FALSE, legend=FALSE, add=TRUE, xpd=xpd)

	if (isTRUE(perimeter)) {
		if (missing(pper) || !is.list(pper)) {
			pper <- list()
		}
		pper$x <- e
		pper$xpd <- xpd
		do.call(lines, pper)
		#lines(e, col=perimeter)
	}

	if (!is.null(box)) {
		if (missing(pbox) || !is.list(pbox)) {
			pbox <- list()
		}
		pbox$x <- box
		pbox$xpd <- xpd
		do.call(lines, pbox)
	}
	invisible(y)
}


setMethod("inset", signature(x="SpatVector"),
	function(x, e, loc="", scale=0.2, background="white", perimeter=TRUE, box=NULL, pper, pbox, ...) {
		.inset(x, e, loc=loc, scale=scale, background=background, perimeter=perimeter, pper=pper, box=box, pbox=pbox, ...)
	}
)


setMethod("inset", signature(x="SpatRaster"),
	function(x, e, loc="", scale=0.2, background="white", perimeter=TRUE, box=NULL, pper, pbox, ...) {
		.inset(x, e, loc=loc, scale=scale, background=background, perimeter=perimeter, pper=pper, box=box, pbox=pbox, ...)
	}
)

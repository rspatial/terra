
RS_locator <- function(n, type, id=FALSE, pch=20, ...) {
# locator that also works in RStudio
# Berry Boessenkool
# https://stackoverflow.com/a/65147220/635245
	on.exit(return(cbind(x, y)))
	x <- y <- NULL
	for (i in seq_len(n)) {
		p <- graphics::locator(1)
		if (is.null(p)) break # ESC
		x <- c(x, p$x)
		y <- c(y, p$y)
		points(x, y, type=type, pch=pch, ...)
		if (id) {
			text(p$x, p$y, labels=i, pos=4, ...)
		}
	}
}

.drawPol <- function(n=1000, id=FALSE, ...) {
	#xy <- graphics::locator(n=1000, type="l", col=col, lwd=lwd, ...)
	#xy <- cbind(xy$x, xy$y)
	xy <- RS_locator(n, "l", id=id, ...)
	xy <- rbind(xy, xy[1,])
	graphics::lines(xy[(length(xy[,1])-1):length(xy[,1]),], ...)
	g <- cbind(1,1,xy,0)
	vect(g, "polygons")
}


.drawLin <- function(n=1000, ...) {
	#xy <- graphics::locator(n=1000, type="l", col=col, lwd=lwd, ...)
	#xy <- cbind(xy$x, xy$y)
	xy <- RS_locator(n, "l", ...)
	g <- cbind(1,1,xy)
	vect(g, "lines")
}


.drawPts <- function(n=1000, ...) {
	#xy <- graphics::locator(n=1000, type="p", col=col, lwd=lwd, ...)
	#xy <- cbind(xy$x, xy$y)
	xy <- RS_locator(n, "p", ...)
	g <- cbind(1:nrow(xy), 1, xy)
	vect(g, "points")
}

.drawExt <- function(...) {
	loc1 <- graphics::locator(n=1, type="p", pch="+", ...)
	loc2 <- graphics::locator(n=1, type="p", pch="+", ...)
	loc <- rbind(unlist(loc1), unlist(loc2))
	e <- c(min(loc[,'x']), max(loc[,'x']), min(loc[,'y']), max(loc[,'y']))
	if (e[1] == e[2]) {
		e[1] <- e[1] - 0.0000001
		e[2] <- e[2] + 0.0000001
	}
	if (e[3] == e[4]) {
		e[3] <- e[3] - 0.0000001
		e[4] <- e[4] + 0.0000001
	}
	p <- rbind(c(e[1], e[3]), c(e[1], e[4]), c(e[2], e[4]), c(e[2], e[3]), c(e[1], e[3]) )
	graphics::lines(p, ...)
	return(ext(e))
}

setMethod("draw", signature(x="character"),
    function(x="extent", col="red", lwd=2, id=FALSE, n=1000, xpd=TRUE, ...){
		x <- match.arg(tolower(x), c("extent", "polygon", "lines", "points"))
		if (x == "extent") {
			.drawExt(col=col, lwd=lwd, xpd=xpd, ...)
		} else if (x == "polygon") {
			.drawPol(n, col=col, lwd=lwd, id=id, xpd=xpd, ...)
		} else if (x == "lines") {
			.drawLin(n, col=col, lwd=lwd, id=id, xpd=xpd, ...)
		} else if (x == "points" || x == "multipoints" ) {
			.drawPts(n, col=col, id=id, xpd=xpd, ...)
		}
	}
)

setMethod("draw", signature(x="missing"),
    function(x="extent", ...){
		draw("extent", ...)
	}
)

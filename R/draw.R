
.drawPol <- function(col="red", lwd=2, ...) {
	xy <- graphics::locator(n=10000, type="l", col=col, lwd=lwd, ...)
	xy <- cbind(xy$x, xy$y)
	xy <- rbind(xy, xy[1,])
	graphics::lines(xy[(length(xy[,1])-1):length(xy[,1]),], col=col, lwd=lwd, ...)
	g <- cbind(1,1,xy,0)
	vect(g, "polygons")
}


.drawLin <- function(col="red", lwd=2, ...) {
	xy <- graphics::locator(n=10000, type="l", col=col, lwd=lwd, ...)
	xy <- cbind(xy$x, xy$y)
	g <- cbind(1,1,xy)
	vect(g, "lines")
}


.drawPts <- function(col="red", lwd=2, ...) {
	xy <- graphics::locator(n=10000, type="p", col=col, lwd=lwd, ...)
	xy <- cbind(xy$x, xy$y)
	g <- cbind(1:nrow(xy), 1, xy)
	vect(g, "points")
}

.drawExt <- function(col="red", lwd=2, ...) {
	loc1 <- graphics::locator(n=1, type="p", pch='+', col=col, ...)
	loc2 <- graphics::locator(n=1, ...)
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
	graphics::lines(p, col=col)
	return(ext(e))
}

setMethod("draw", signature(x="character"),
    function(x="extent", col="red", lwd=2, ...){ 
		objtypes <- c("extent", "polygon", "line", "points")
		i <- pmatch(tolower(x), objtypes)
		if (is.na(i)) {
			stop("invalid object type")
		} else if (i < 1) {
			stop("ambiguous object type")
		}
		x <- objtypes[i]
		if (x == "extent") {
			.drawExt(col, lwd, ...)
		} else if (x == "polygon") {
			.drawPol(col, lwd, ...)
		} else if (x == "lines") {
			.drawLin(col, lwd, ...)
		} else if (x == "points") {
			.drawPts(col, lwd, ...)
		} 
	}
)

setMethod("draw", signature(x="missing"),
    function(x="extent", col="red", lwd=2, ...){ 
		draw("extent", col, lwd, ...)
	}
)

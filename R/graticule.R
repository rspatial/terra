# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date : November 2022
# Version 1
# License GPL v3

graticule <- function(lon=30, lat=30, crs="") {

	if (crs == "") {
		crs <- "+proj=longlat"
	}
	
	if (length(lon) == 1) {
		lon <- seq(-180, 180, lon)
	} else {
		lon <- sort(lon)
	}
	
	if (length(lat) == 1) {
		lat <- seq(-90, 90, lat)
	} else {
		lat <- sort(lat)
	}

	interval <- 100 #should depend on extent
	
	y <- cbind(rep(1:length(lon), each=2), rep(lon, each=2), range(lat))
	vy <- vect(y, "lines", crs="+proj=longlat")
	vy <- densify(vy, interval*1000, TRUE)
	values(vy) <- data.frame(h=FALSE, lon=lon)

	rlon <- range(lon)
	rlon <- c(rlon[1], rlon[1] + (rlon[2]-rlon[1])/2, rlon[2])
	x <- cbind(rep(1:length(lat), each=3), rlon, rep(lat, each=3))
	vx <- vect(x, "lines", crs="local")
	vx <- densify(vx, interval/110, TRUE)
	values(vx) <- data.frame(h=TRUE, lat=lat)
	crs(vx, warn=FALSE) <- "+proj=longlat"
	
	v <- rbind(vy, vx)
	
	e <- as.polygons(ext(v), crs=crs(v))
	e <- densify(e, interval*1000, TRUE)
	
	if (!is.lonlat(crs, FALSE, FALSE)) {
		v <- project(v, crs)
		e <- project(e, crs)
	}

	g <- new("SpatGraticule")
	g@pntr <- v@pntr
	g@box <- e@pntr 
	g
}


#setMethod("project", signature(x="SpatGraticule"),
#	function(x, y) {
#		v <- vect(x)
#		v <- project(v, y)
#		x@pntr <- v@pntr
#		v@pntr <- x@box
#		v <- project(v, y)
#		x@box <- v@pntr
#		x
#	}
#)

setMethod("crop", signature(x="SpatGraticule"),
	function(x, y) {
		v <- vect(x)

		v <- crop(v, y)
		x@pntr <- v@pntr

		v@pntr <- x@box
		v <- crop(v, y)
		x@box <- v@pntr
		x
	}
)

setMethod("erase", signature(x="SpatGraticule", y="SpatVector"),
	function(x, y) {

		v <- vect(x)
		y <- project(y, v)

		v <- erase(v, y)
		x@pntr <- v@pntr

		v@pntr <- x@box
		v <- erase(v, y)
		x@box <- v@pntr
		x
	}
)


grat_tick <- function(x, direction, tick.length) {

	y <- crds(x)
	if (direction == "left") {
		y[,1] <- y[,1] - tick.length 
	} else if (direction == "right") {
		y[,1] <- y[,1] + tick.length 
	} else if (direction == "top") {
		y[,2] <- y[,2] + tick.length 
	} else {
		y[,2] <- y[,2] - tick.length 
	}
	lines(x, vect(y), xpd=TRUE)

}


grat_labels <- function(x, retro, atlon, atlat, labloc, cex, col, offlon, offlat, vfont=NULL, font=NULL, tickmarks=FALSE, ...) {

	v <- vect(x)
	if (tickmarks) {
		cw <- graphics::strwidth("A")
		offlat <- 2 * offlat
		offlon <- 2 * offlon
	}
	
	left <- right <- top <- bottom <- FALSE
	labloc <- rep_len(labloc, 2)
	if (!is.na(labloc[1])) {
		if (labloc[1] == 1) {
			bottom <- TRUE
		} else if (labloc[1] == 2) {
			top <- TRUE
		} else {
			bottom <- TRUE
			top <- TRUE	
		}
	}
	if (!is.na(labloc[2])) {
		if (labloc[2] == 1) {
			left <- TRUE
		} else if (labloc[2] == 2) {
			right <- TRUE
		} else {
			left <- TRUE
			right <- TRUE	
		}
	}
	if (left || right) {
		s <- sign(offlat)+2
		x <- v[v$h, ]
		g <- geom(x)
		if (retro) {
			labs <- retro_labels(x$lat, lat=TRUE)
		} else {
			labs <- paste0(x$lat,"\u00B0")		
		}	
		if (!is.null(atlat)) {
			atlat <- round(stats::na.omit(atlat))
			atlat <- atlat[(atlat > 0) & (atlat <= length(labs))] 
		}
		if (left) {
			a <- vect(g[match(unique(g[,1]), g[,1]), 3:4])
			a$labs <- labs
			if (!is.null(atlat)) {
				a <- a[atlat, ]
			}
			pos <- c(4,0,2)[s]
			if (pos == 0) pos <- NULL
			text(a, labels="labs", pos=pos, offset=abs(offlat), cex=cex, halo=TRUE, xpd=TRUE, col=col)
			if (tickmarks) {
				grat_tick(a, "left", cw * abs(offlat))
			}
		}
		if (right) {
			g <- g[nrow(g):1, ]
			a <- vect(g[match(unique(g[,1]), g[,1]), 3:4])
			a$labs <- rev(labs)
			if (!is.null(atlat)) {
				a <- a[atlat, ]
			}
			pos <- c(2,0,4)[s]			
			if (pos == 0) pos <- NULL
			text(a, labels="labs", pos=pos, offset=abs(offlat), cex=cex, halo=TRUE, xpd=TRUE, col=col)
			if (tickmarks) {
				grat_tick(a, "right", cw * abs(offlat))
			}
		}
	}
	if (top || bottom) {
		s <- sign(offlon)+2
		x <- v[!v$h, ]
		g <- geom(x)
		if (retro) {
			labs <- retro_labels(x$lon, lat=FALSE)
		} else {
			labs <- paste0(x$lon,"\u00B0")		
		}
		if (!is.null(atlon)) {
			atlon <- round(stats::na.omit(atlon))
			atlon <- atlon[(atlon > 0) & (atlon <= length(labs))] 
		}
		if (bottom) {
			a <- vect(g[match(unique(g[,1]), g[,1]), 3:4])
			a$labs <- labs
			if (!is.null(atlon)) {
				a <- a[atlon, ]
			}
			pos <- c(3,0,1)[s]
			if (pos == 0) pos <- NULL
			text(a, labels="labs", pos=pos, offset=abs(offlon), cex=cex, halo=TRUE, xpd=TRUE, col=col)
			if (tickmarks) {
				grat_tick(a, "bottom", cw * abs(offlon))
			}
		}
		if (top) {
			g <- g[nrow(g):1, ]
			a <- vect(g[match(unique(g[,1]), g[,1]), 3:4])
			a$labs <- labs
			if (!is.null(atlon)) {
				a <- a[atlon, ]
			}
			pos <- c(1,0,3)[s]
			if (pos == 0) pos <- NULL
			text(a, labels="labs", pos=pos, offset=abs(offlon), cex=cex, halo=TRUE, xpd=TRUE, col=col)
			if (tickmarks) {
				grat_tick(a, "top", cw * abs(offlon))
			}
		}
	}	
}	



setMethod("plot", signature(x="SpatGraticule", y="missing"),
	function(x, y, background=NULL, col="black", mar=NULL, labels=TRUE, retro=FALSE, lab.loc=c(1,1), lab.lon=NULL, lab.lat=NULL, lab.cex=.65, lab.col="black", off.lat=0.25, off.lon=0.25, box=FALSE, box.col="black", tickmarks=FALSE, add=FALSE, ...) {
		b <- vect()
		b@pntr <- x@box
		if (!is.null(mar)) mar <- rep_len(mar, 4)
		if (!is.null(background)) {
			plot(b, col=background, border=NA, axes=FALSE, mar=mar, add=add)
		} else {
			if (!add) {
				plot(ext(b), border=NA, axes=FALSE, mar=mar)
			}
		}
		lines(x, col=col, ...)
		if (box) {
			lwd <- list(...)$lwd
			if (is.null(lwd)) lwd=1;
			lines(b, lty=1, col=box.col, lwd=lwd)
		}
		if (labels) grat_labels(x, isTRUE(retro[1]), lab.lon, lab.lat, lab.loc, lab.cex[1], lab.col[1], off.lon, off.lat, tickmarks=tickmarks, ...)
	}
)


setMethod("lines", signature(x="SpatGraticule"),
	function(x, ...) {
		lines(vect(x), ...)
	}
)



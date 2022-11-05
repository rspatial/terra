# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date : November 2022
# Version 1
# License GPL v3

graticule <- function(lat, lon=lat, crs="+proj=longlat") {

	interval <- 100

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
	crs(vx) <- "+proj=longlat"

	v <- rbind(vy, vx)
	e <- as.polygons(ext(v), crs=crs(v))
	e <- densify(e, interval*1000, TRUE)
	
	if (crs != "+proj=longlat") {
		v <- project(v, crs)
		e <- project(e, crs)
	}

	g <- new("SpatGraticule")
	g@ptr <- v@ptr 
	g@box <- e@ptr 
	g
}


setMethod("project", signature(x="SpatGraticule"),
	function(x, y) {
		v <- vect()
		v@ptr <- x@ptr
		v <- project(v, y)
		x@ptr <- v@ptr
		v@ptr <- x@box
		v <- project(v, y)
		x@box <- v@ptr
		x
	}
)

setMethod("crop", signature(x="SpatGraticule"),
	function(x, y) {
		v <- vect()

		v@ptr <- x@ptr
		v <- crop(v, y)
		x@ptr <- v@ptr

		v@ptr <- x@box
		v <- crop(v, y)
		x@box <- v@ptr
		x
	}
)


setMethod("plot", signature(x="SpatGraticule", y="missing"),
	function(x, y, background=NULL, col="black", mar=c(1,1,1,1), ...) {
		v <- vect()
		if (!is.null(background)) {
			v@ptr <- x@box
			plot(v, col=background, border=NA, axes=FALSE, mar=mar)
			v@ptr <- x@ptr
		} else {
			v@ptr <- x@ptr
			plot(v, axes=FALSE, type="n", mar=mar)
		}
		lines(v, col=col, ...)
	}
)

setMethod("lines", signature(x="SpatGraticule"),
	function(x, ...) {
		v <- vect()
		v@ptr <- x@ptr
		lines(v, ...)
	}
)



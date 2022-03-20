# Author: Robert J. Hijmans 
# Date : October 2018
# Version 1.0
# License GPL v3

## from stars
#stars:::st_as_raster is used
#setAs("stars", "SpatRaster") is provided by stars via st_as_raster


from_stars <- function(from) {

	isProxy <- inherits(from, "stars_proxy")
	natts <- length(from)
	#from[i] recursion does not work with proxy
	if (!isProxy && (natts > 1)) { # not sure what attributes represent
		ra <- list()
		for (i in 1:natts) {
			ra[[i]] <- from_stars(from[i])
		}
		if (all(sapply(ra, function(i) inherits(i, "SpatRaster")))) {
			nl <- sapply(ra, nlyr)
			ra <- rast(ra)
			nms <- names(ra)
			names(ra) <- paste(rep(names(from), nl), nms, sep="_")
		} else 	if (all(sapply(ra, function(i) inherits(i, "SpatRasterDataset")))) {
			ra <- do.call(c, ra)
		} else {
			ra <- lapply(ra, function(i) if (!inherits(i, "SpatRasterDataset")) {sds(i)} else {i})
			ra <- do.call(c, ra)
		}
		return(ra)
	}

	dims <- attr(from, "dimensions")
	dd <- dim(from)

	# x, y
	hasBands <- "band" %in% names(dd)
	hasTime <- "time" %in% names(dd)
	timev <- NULL
	if (hasTime) {
		tim <- dims$time$offset
		tseq <- dims$time$from:dims$time$to
		if (dims$time$refsys == "Date") {
			timev <- as.Date(tim) + tseq
		} else { # for now
			timev <- tseq
		}
	}

	# no time or variables
	if (length(dd) - hasBands == 2) {
		return( methods::as(from, "SpatRaster"))
	}


	# time, perhaps bands or variables
	if (length(dd) - (hasTime + hasBands) == 2) {
		r <- methods::as(from, "SpatRaster")
		if (hasBands) {
			timev <- rep(timev, each=dd["band"])
		} 
		time(r) <- timev
		return(r)
	}

	if (isProxy) {
		# currently not setting time dim here
		if (natts > 1) {
			ff <- sapply(from, function(i) from[i][[1]])
			s <- sds(ff)
			names(s) <- names(from) 
		} else {
			f <- from[[1]]
			s <- sds(f)
			nms <- names(dd)[3+hasBands]
			if (!is.na(nms)) {
				names(s) <- paste(nms, 1:length(s), sep="-")
			}
		}
		return(s)
	}

	xmin <- dims$x$offset
	nc <- dims$x$to
	xmax <- xmin + nc * dims$x$delta
	ymax <- dims$y$offset
	nr <- dims$y$to
	ymin <- ymax + nr * dims$y$delta

	from <- from[[1]]
	rr <- list()
	if (hasTime && hasBands) {
		for (i in 1:dd[5]) {
			x <- from[,,,,i]
			r <- rast(ncols=nc, nrows=nr, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, crs=dims$x$refsys$wkt, nlyr=dd["band"] * dd["time"])
			time(r) <- rep(timev, each=dd["band"])
			bandnames <- rep(paste("band", 1:dd["band"], sep="-"), length(timev))
			names(r) <- paste(bandnames, rep(timev, each=dd["band"]), sep="_")
			r <- setValues(r, as.vector(x))
			rr[[i]] <- r
		}
	} else { #if (hasTime || hasBands) {
		for (i in 1:dd[4]) {
			x <- from[,,,i]
			r <- rast(ncols=nc, nrows=nr, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, crs=dims$x$refsys$wkt, nlyr=dim(x)[3], time=timev)
			if (hasBands) {
				names(r) <- paste("band", 1:dd["band"], sep="-")
			} else {
				names(r) <- timev
			}
			rr[[i]] <- setValues(r, x)
		}
	} 
	s <- sds(rr)
	names(s) <- paste(names(dd)[4], 1:length(s), sep="-")
	s
}



setAs("stars", "SpatRasterDataset",
	function(from) {
		from_stars(from) 
	}
)

setAs("ggmap", "SpatRaster", 
	function(from) {
		b <- attr(from, "bb")
		e <- ext(b$ll.lon, b$ur.lon, b$ll.lat, b$ur.lat)
		r <- rast(nrows=nrow(from), ncols=ncol(from), ext=e, nlyr=3, crs="epsg:4326", names=c("red", "green", "blue"))
		values(r) <- t(grDevices::col2rgb(from))
		RGB(r) <- 1:3
		r
	}
)



as.list.SpatRaster <- function(x, ...) {
	lapply(1:nlyr(x), function(i) x[[i]])
}
setMethod("as.list", signature(x="SpatRaster"), as.list.SpatRaster)


 
# create a "grDevices::raster" (small r) object for use with the rasterImage function
# NOT a raster::Raster* object
setMethod("as.raster", signature(x="SpatRaster"), 
	function(x, maxcell=500000, col) {
		if (missing(col)) {
			col <- rev(grDevices::terrain.colors(255))
		}
		x <- spatSample(x, maxcell, method="regular", as.raster=TRUE)
		x <- as.matrix(x, wide=TRUE)
		r <- range(x, na.rm=TRUE)
		x <- (x - r[1])/ (r[2] - r[1])
		x <- round(x * (length(col)-1) + 1)
		x[] <- col[x]
		as.raster(x)
	} 
)



.as.image <- function(x, maxcells=10000) {
	x <- spatSample(x, size=maxcells, method="regular", as.raster=TRUE)
	X <- xFromCol(x, 1:ncol(x))
	Y <- yFromRow(x, nrow(x):1)
	Z <- t(as.matrix(x, wide=TRUE)[nrow(x):1,]) 
	list(x=X, y=Y, z=Z)
}
 
 
setMethod("as.polygons", signature(x="SpatRaster"), 
	function(x, trunc=TRUE, dissolve=TRUE, values=TRUE, na.rm=TRUE, extent=FALSE) {
		p <- methods::new("SpatVector")
		if (extent) {
			p@ptr <- x@ptr$dense_extent(FALSE, FALSE)
			x <- messages(x, "as.polygons")
		} else {
			opt <- spatOptions()
			p@ptr <- x@ptr$as_polygons(trunc[1], dissolve[1], values[1], na.rm[1], opt)
			x <- messages(x, "as.polygons")
			if (values) {
				ff <- is.factor(x)
				if (dissolve) {
					ff <- ff[[1]]
				}
				if (any(ff)) {
					ff <- which(ff)
					cgs <- cats(x)
					for (f in ff) {
						cg <- cgs[[f]]
						i <- match(unlist(p[[f]]), cg[,1])
						act <- activeCat(x, f)
						p[[f]] <- cg[i, act+1]
					}
				}
			}
		}
		messages(p, "as.polygons")
	}
)

setMethod("as.lines", signature(x="SpatRaster"), 
	function(x) {
		p <- methods::new("SpatVector")
		opt <- spatOptions()
		p@ptr <- x@ptr$as_lines(opt)
		messages(p, "as.lines")
	}
)


setMethod("as.polygons", signature(x="SpatExtent"), 
	function(x, crs="") {
		p <- methods::new("SpatVector")
		crs <- character_crs(crs, "as.polygons")
		p@ptr <- SpatVector$new(x@ptr, crs)
		messages(p, "as.polygons")
	}
)

setMethod("as.lines", signature(x="SpatExtent"), 
	function(x, crs="") {
		crs <- character_crs(crs, "lines")
		as.lines(as.polygons(x, crs))
	}
)


setMethod("as.points", signature(x="SpatExtent"), 
	function(x, crs="") {
		#vect(do.call(cbind, x@ptr$as.points()), "points", crs=crs)
		as.points(as.polygons(x, crs))
	}
)


setMethod("as.lines", signature(x="SpatVector"), 
	function(x) {
		x@ptr <- x@ptr$as_lines()
		messages(x, "as.lines")
	}
)

setMethod("as.polygons", signature(x="SpatVector"), 
	function(x) {
		x@ptr <- x@ptr$polygonize()
		messages(x, "as.polygons")
	}
)

setMethod("as.points", signature(x="SpatVector"), 
	function(x, multi=FALSE, skiplast=TRUE) {
		opt <- spatOptions()
		x@ptr <- x@ptr$as_points(multi, skiplast)
		messages(x, "as.points")
	}
)


setMethod("as.points", signature(x="SpatRaster"), 
	function(x, values=TRUE, na.rm=TRUE) {
		p <- methods::new("SpatVector")
		opt <- spatOptions()
		p@ptr <- x@ptr$as_points(values, na.rm, opt)
		x <- messages(x, "as.points")

		if (values) {
			ff <- is.factor(x)
			if (any(ff)) {
				ff <- which(ff)
				levs <- levels(x)
				for (f in ff) {
					facts <- levs[[f]]
					v <- factor(unlist(p[[f]], use.names=FALSE), levels=(1:length(facts))-1)
					levels(v) <- facts
					p[[f]] <- as.character(v)
				}
			}
		}
		messages(p, "as.points")
	}
)

# mode argument is ignored as mode=mode gave an error on R-devel
setMethod("as.vector", signature(x="SpatExtent"), 
	function(x, mode="any") {
		v <- x@ptr$vector
		names(v) <- c("xmin", "xmax", "ymin", "ymax")
		if (mode == "list") {
			v <- as.list(v)
		}
		v
	}
)


setMethod("as.character", signature(x="SpatExtent"), 
	function(x) {
		e <- as.vector(x)
		paste0("ext(", paste(e, collapse=", "), ")")
	}
)

setMethod("as.vector", signature(x="SpatRaster"), 
	function(x, mode="any") {
		values(x, FALSE)
	}
)

as.matrix.SpatRaster <- function(x, ...) {
	if (!hasValues(x)) {
		error("as.matrix", "SpatRaster has no cell values")
	}
	wide <- isTRUE(list(...)$wide)
	if (wide) {
		if (nlyr(x) > 1) {
			m <- values(x, mat=TRUE)
			m <- lapply(1:ncol(m), function(i) {
				matrix(m[,i], nrow=nrow(x), byrow=TRUE)
				})
			m <- do.call(cbind, m)
		} else {
			m <- matrix(values(x, mat=FALSE),nrow=nrow(x),byrow=TRUE)
		}
	} else {
		m <- values(x, mat=TRUE)
	}
	m
}
setMethod("as.matrix", signature(x="SpatRaster"), as.matrix.SpatRaster)


as.data.frame.SpatRaster <- function(x, row.names=NULL, optional=FALSE, xy=FALSE, cells=FALSE, na.rm=TRUE, ...) {
#	dots <- list(...) 
#	xy <- isTRUE(dots$xy) 
#	cells <- isTRUE(dots$cells)
#	na.rm <- isTRUE(dots$na.rm)

	d <- NULL
	if (xy) {
		d <- xyFromCell(x, 1:ncell(x))
	} 
	if (cells) {
		d <- cbind(cell=1:ncell(x), d)
	}
	if (is.null(d)) {
		d <- values(x, dataframe=TRUE, ... )
	} else {
		d <- data.frame(d)
		d <- cbind(d, values(x, dataframe=TRUE), ...)
	}
	if (na.rm) {
		d <- stats::na.omit(d) 
		attr(d, "na.action") <- NULL
	}
	d
}
setMethod("as.data.frame", signature(x="SpatRaster"), as.data.frame.SpatRaster)



setAs("SpatRaster", "data.frame", 
	function(from) {
		as.data.frame(from)
	}
)

setAs("SpatVector", "data.frame", 
	function(from) {
		as.data.frame(from)
	}
)

setMethod("as.array", signature(x="SpatRaster"), 
	function(x) {
		dm <- dim(x)
		x <- values(x, TRUE)
		a <- array(NA, dm)
		for (i in 1:dm[3]) {
			a[,,i] <- matrix(x[,i], nrow=dm[1], byrow=TRUE)
		}
		a
	}
)


# to sf from SpatVector
# available in sf
#.v2sf <- function(from) {
#	txt <- 'sf::st_as_sf(as.data.frame(from, geom=TRUE), wkt="geometry", crs=from@ptr$get_crs("wkt"))'
#	eval(parse(text = txt))
#}

# sf bbox
.ext_from_sf <- function(from) {
	sfi <- attr(from, "sf_column")
	geom <- from[[sfi]]
	e <- attr(geom, "bbox")
	ext(e[c(1,3,2,4)])
}


.from_sf <- function(from) {
	sfi <- attr(from, "sf_column")
	geom <- from[[sfi]]
	crs <- attr(geom, "crs")$wkt
	if (is.na(crs)) crs <- ""
	#geom <- st_as_text(geom)
	#v <- vect(geom, crs=crs)
	v <- vect()
	v@ptr <- v@ptr$from_hex(sf::rawToHex(sf::st_as_binary(geom)), crs)
	if (ncol(from) > 1) {
		from[[sfi]] <- NULL
		values(v) <- as.data.frame(from)
	}
	v
}


setAs("sf", "SpatVector", 
	function(from) {
		v <- try(.from_sf(from), silent=TRUE)
		if (inherits(v, "try-error")) {
			error("as,sf", "coercion failed. You can try coercing via a Spatial* (sp) class")
		} 
		v
	}
)

.from_sfc <- function(from) {
	v <- vect()
	v@ptr <- v@ptr$from_hex(sf::rawToHex(sf::st_as_binary(from)), "")
	crs(v) <- attr(from, "crs")$wkt
	v
}


setAs("sfc", "SpatVector", 
	function(from) {
		v <- try(.from_sfc(from), silent=TRUE)
		if (inherits(v, "try-error")) {
			error("as,sfc", "coercion failed. You can try coercing via a Spatial* (sp) class")
		} 
		v
	}
)


setAs("sfg", "SpatVector", 
	function(from) {
		v <- try(.from_sfc(from), silent=TRUE)
		if (inherits(v, "try-error")) {
			error("as,sfg", "coercion failed. You can try coercing via a Spatial* (sp) class")
		}
		v
	}
)

setAs("XY", "SpatVector", 
	function(from) {
		v <- try(.from_sfc(from), silent=TRUE)
		if (inherits(v, "try-error")) {
			error("as,sfc", "coercion failed. You can try coercing via a Spatial* (sp) class")
		} 
		v
	}
)

setAs("im", "SpatRaster", 
	function(from) {
		r <- rast(nrows=from$dim[1], ncols=from$dim[2], xmin=from$xrange[1], xmax=from$xrange[2], ymin=from$yrange[1], ymax=from$yrange[2], crs="")
		values(r) <- from$v
		flip(r, direction="vertical")
	}
)


setAs("SpatVector", "Spatial", 
	function(from) {
		g <- geom(from, df=TRUE)
		geom(g, values(from), geomtype(from), crs(from))
	}
)


setAs("Spatial", "SpatVector", 
	function(from) {
		g <- geom(from, df=TRUE)
		colnames(g)[1] <- "id"
		if (inherits(from, "SpatialPolygons")) {
			vtype <- "polygons"
			if ("cump" %in% colnames(g)) {
				g <- g[,c(1,2,5,6,4)]
			}
		} else if (inherits(from, "SpatialLines")) {
			vtype <- "lines"
			if ("cump" %in% colnames(g)) {
				g <- g[,colnames(g) != "cump"]
			}
		} else {
			vtype <- "points"
			g <- cbind(g[,1,drop=FALSE], part=1:nrow(g), g[,2:3,drop=FALSE])
		}
		if (methods::.hasSlot(from, "data")) {
			v <- vect(g, vtype, from@data, crs(from))
		} else {
			v <- vect(g, vtype, crs=crs(from))
		}
		return(v)
	}
)




setAs("SpatialGrid", "SpatRaster", 
	function(from){
		prj <- attr(from@proj4string, "comment")
		if (is.null(prj)) prj <- from@proj4string@projargs
		b <- rast(ext=as.vector(t(from@bbox)), crs=prj)
		if (inherits(from, "SpatialGridDataFrame")) {
			dim(b) <- c(from@grid@cells.dim[2], from@grid@cells.dim[1], ncol(from@data))		
			names(b) <- colnames(from@data)
			b <- setValues(b, as.matrix(from@data))
		} else {
			dim(b) <- c(from@grid@cells.dim[2], from@grid@cells.dim[1])		
		}
		b
	}
)

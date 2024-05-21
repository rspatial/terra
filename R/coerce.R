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



as.list.SpatRaster <- function(x, geom=NULL, ...) {
	if (!is.null(geom)) {
		e <- as.vector(ext(x))
		d <- crs(x, describe=TRUE)
		if (!(is.na(d$authority) || is.na(d$code))) {
			crs <- paste0(d$authority, ":", d$code)
		} else {
			crs <- gsub("\n[ ]+", "", crs(x))
		}
		list(
			ncols=ncol(x),
			nrows=nrow(x),
			nlyrs=nlyr(x),
			xmin=e[1],
			xmax=e[2],
			ymin=e[3],
			ymax=e[4],
			xres=xres(x),
			yres=yres(x),
			nms = paste(names(x), collapse="', '"),
			units = paste(units(x), collapse="', '"),
			time = paste(time(x), collapse="', '"),
			crs=crs
		)
	} else {
		lapply(1:nlyr(x), function(i) x[[i]])
	}
}
setMethod("as.list", signature(x="SpatRaster"), as.list.SpatRaster)


as.list.SpatRasterCollection <- function(x, ...) {
	out <- lapply(1:length(x), function(i) x[i])
	names(out) <- names(x)
	out
}
setMethod("as.list", signature(x="SpatRasterCollection"), as.list.SpatRasterCollection)


as.list.SpatRasterDataset <- function(x, ...) {
	out <- lapply(1:length(x), function(i) x[i])
	names(out) <- names(x)
	out
}
setMethod("as.list", signature(x="SpatRasterDataset"), as.list.SpatRasterDataset)


as.list.SpatVectorCollection <- function(x, ...) {
	out <- lapply(1:length(x), function(i) x[i])
	names(out) <- names(x)
	out
}
setMethod("as.list", signature(x="SpatVectorCollection"), as.list.SpatVectorCollection)


# create a "grDevices::raster" (small r) object for use with the rasterImage function
# NOT a raster::Raster* object
setMethod("as.raster", signature(x="SpatRaster"),
	function(x, maxcell=500000, col) {
		if (missing(col)) {
			#col <- rev(grDevices::terrain.colors(255))
			col <- .default.pal()
		}
		x <- spatSample(x, maxcell, method="regular", as.raster=TRUE, warn=FALSE)
		x <- as.matrix(x, wide=TRUE)
		r <- range(x, na.rm=TRUE)
		x <- (x - r[1])/ (r[2] - r[1])
		x <- round(x * (length(col)-1) + 1)
		x[] <- col[x]
		as.raster(x)
	}
)



.as.image <- function(x, maxcells=10000) {
	x <- spatSample(x, size=maxcells, method="regular", as.raster=TRUE, warn=FALSE)
	X <- xFromCol(x, 1:ncol(x))
	Y <- yFromRow(x, nrow(x):1)
	Z <- t(as.matrix(x, wide=TRUE)[nrow(x):1,])
	list(x=X, y=Y, z=Z)
}


get_labels <- function(x, p, dissolve=FALSE) {
	ff <- is.factor(x)
	if (dissolve) {
		ff <- ff[[1]]
	}
	if (any(ff)) {
		ff <- which(ff)
		cgs <- levels(x)
		for (f in ff) {
			cg <- cgs[[f]]
			i <- match(unlist(p[[f]]), cg[,1])
			p[[f]] <- cg[i, 2]
		}
	}
	p
}


setMethod("as.polygons", signature(x="SpatRaster"),
	function(x, round=TRUE, aggregate=TRUE, values=TRUE, na.rm=TRUE, na.all=FALSE, extent=FALSE, digits=0, ...) {
		if (isFALSE(list(...)$dissolve)) {
			aggregate <- FALSE
		}
	
		p <- methods::new("SpatVector")
		if (extent) {
			p@ptr <- x@ptr$dense_extent(FALSE, FALSE)
			x <- messages(x, "as.polygons")
		} else {
			opt <- spatOptions()
			p@ptr <- x@ptr$as_polygons(round[1], aggregate[1], values[1], na.rm[1], na.all[1], digits, opt)
			x <- messages(x, "as.polygons")
			if (values) {
				p <- get_labels(x, p, aggregate[1])
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
		crs <- character_crs(crs, "as.lines")
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
	function(x, extent=FALSE) {
		if (extent) {
			as.polygons(ext(x), crs=crs(x))
		} else {
			x@ptr <- x@ptr$polygonize()
			messages(x, "as.polygons")
		}
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
	function(x, values=TRUE, na.rm=TRUE, na.all=FALSE) {
		p <- methods::new("SpatVector")
		opt <- spatOptions()
		p@ptr <- x@ptr$as_points(values, na.rm, na.all, opt)
		x <- messages(x, "as.points")

		if (values) {
			p <- get_labels(x, p, FALSE)
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

as.matrix.SpatRaster <- function(x, wide=FALSE, ...) {
	if (!hasValues(x)) {
		error("as.matrix", "SpatRaster has no cell values")
	}
#	wide <- isTRUE(list(...)$wide)
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




as.matrix.SpatExtent <- function(x, ...) {
	v <- matrix(as.vector(x), ncol=2, byrow=TRUE)
	colnames(v) <- c("min", "max")
	v
}
setMethod("as.matrix", signature(x="SpatExtent"), as.matrix.SpatExtent)



as.data.frame.SpatRaster <- function(x, row.names=NULL, optional=FALSE, xy=FALSE, cells=FALSE, time=FALSE,na.rm=NA, wide=TRUE, ...) {

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
	if (is.na(na.rm)) {
		cols <- (1 + cells + xy * 2):ncol(d)
		i <- rowSums(is.na(d[,cols,drop=FALSE])) < length(cols)
		d <- d[i,,drop=FALSE]
	} else if (isTRUE(na.rm)) {
		d <- stats::na.omit(d)
		attr(d, "na.action") <- NULL
	}
	if (!wide) {
		nr <- nrow(d)
		if (!(xy || cells)) {
			d <- data.frame(layer=rep(names(x), each=nr), values=as.vector(as.matrix(d)))
		} else {
			idv <- NULL
			if (xy) idv <- c("x", "y", idv)
			if (cells) idv <- c("cell", idv)
			add <- d[idv]
			rownames(add) <- NULL
			for (i in 1:length(idv)) {
				d[[idv[i]]] <- NULL
			}
			d <- data.frame(add, layer=rep(names(x), each=nr), values=as.vector(as.matrix(d)))
			nms <- names(x)
#			d <- stats::reshape(d, direction="long", idvar=idv, varying=nms, v.names="values")
#			d$time <- nms[d$time]
#			names(d)[names(d) == "time"] <- "layer"
		}
		rownames(d) <- NULL
		if (time) {
		#	d$time <- NULL
		#	vals <- d$values
		#	d$values <- NULL
			d$time <- rep(time(x), each=nr)
		#	d$values <- vals
		}
	} else if (time && has.time(x)) {
		tm <- as.character(time(x))
		nc <- ncol(d)
		colnames(d)[(1+nc-length(tm)):nc] <- tm
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


.from_sf <- function(from, geom, sfi) {
	crs <- attr(geom, "crs")$wkt
	if (is.na(crs)) crs <- ""
	#geom <- st_as_text(geom)
	#v <- vect(geom, crs=crs)
	v <- vect()
	v@ptr <- v@ptr$from_hex(sf::rawToHex(sf::st_as_binary(geom)), crs)
	v <- messages(v, "SpatVector from sf")
	if (ncol(from) > 1) {
		from[[sfi]] <- NULL
		values(v) <- as.data.frame(from)
	}
	v
}

.svc_from_sf <- function(from) {
	sfi <- attr(from, "sf_column")
	geom <- from[[sfi]]
	crs <- attr(geom, "crs")$wkt
	if (is.na(crs)) crs <- ""
	#geom <- st_as_text(geom)
	#v <- vect(geom, crs=crs)
	v <- svc()
	v@ptr <- v@ptr$from_hex_col(sf::rawToHex(sf::st_as_binary(geom)), crs)
	#if (ncol(from) > 1) {
	#	from[[sfi]] <- NULL
	#	values(v) <- as.data.frame(from)
	#}
	v
}

setAs("sf", "SpatRaster",
	function(from) {
		e <- ext(from)
		rast(e)
	}
)


setAs("sf", "SpatVector",
	function(from) {
		sfi <- attr(from, "sf_column")
		if (is.null("sfi")) {
			error("as,sf", "the object has no sf_column")
		}
		geom <- from[[sfi]]
		if (inherits(geom, "list")) {
			error("as,sf", "the geometry column is not valid (perhaps first load the sf package)")
		}
		v <- try(.from_sf(from, geom, sfi), silent=FALSE)
		if (inherits(v, "try-error")) {
			error("as,sf", "coercion failed. You can try coercing via a Spatial* (sp) class")
		}
		v
	}
)

.from_sfc <- function(from) {
	v <- vect()
	v@ptr <- v@ptr$from_hex(sf::rawToHex(sf::st_as_binary(from)), "")
	crs(v, warn=FALSE) <- attr(from, "crs")$wkt
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
		units(r) <- from$units$singular
		if (from$units$multiplier != 1) {
			r <- r * from$units$multiplier
		}
		flip(r, direction="vertical")
	}
)


setAs("SpatVector", "Spatial",
	function(from) {
		hasmethod <- suppressWarnings("geom,data.frame-method" %in% utils::methods("geom"))
		if (!hasmethod) {
			error("coerce", "first run 'library(raster)' to coerce a SpatVector to a Spatial object" )
		}
		g <- geom(from, df=TRUE)
		geom(g, values(from), geomtype(from), as.character(crs(from, proj=TRUE)))
	}
)


geom_SpatialPolygons <- function(x) {
	nobs <- length(x@polygons)
	objlist <- vector(mode = "list", length = nobs)
	for (i in 1:nobs) {
		nsubobs <- length(x@polygons[[i]]@Polygons)
		ps <- list()
		last <- 0
		for (j in 1:nsubobs) {
			if (!x@polygons[[i]]@Polygons[[j]]@hole) {
				last <- last + 1
				hole <- 0
			} else {
				hole <- max(1, last)
			}
			ps[[j]] <- cbind(j, x@polygons[[i]]@Polygons[[j]]@coords, hole)
		}
		objlist[[i]] <- cbind(i, do.call(rbind, ps))
	}
	do.call(rbind, objlist)
}

geom_SpatialLines <- function(x) {
	nobs <- length(x@lines)
	objlist <- vector(mode = "list", length = nobs)
	for (i in 1:nobs) {
		nsubobj <- length(x@lines[[i]]@Lines)
		ps <- lapply(1:nsubobj, function(j) cbind(j, x@lines[[i]]@Lines[[j]]@coords))
		objlist[[i]] <- cbind(i, do.call(rbind, ps))
	}
	do.call(rbind, objlist)
}


setAs("Spatial", "SpatVector",
	function(from) {
		if (inherits(from, "SpatialPolygons")) {
			g <- geom_SpatialPolygons(from)
			vtype <- "polygons"
		} else if (inherits(from, "SpatialLines")) {
			g <- geom_SpatialLines(from)
			vtype <- "lines"
		} else if (inherits(from, "SpatialPoints")) {
			g <- from@coords[,1:2,drop=FALSE]
			vtype <- "points"
		} else {
			error("coerce", "cannot coerce this object to a SpatVector")
		}
		#the below can change the proj-string when going back to sp
		#crs <- attr(from@proj4string, "comment")
		#if (is.null(crs)) 
		crs <- from@proj4string@projargs
		if (methods::.hasSlot(from, "data")) {
			vect(g, vtype, from@data, crs=crs)
		} else {
			vect(g, vtype, crs=crs)
		}
	}
)


setAs("SpatialGrid", "SpatRaster",
	function(from){
		prj <- attr(from@proj4string, "comment")
		if (is.null(prj)) prj <- from@proj4string@projargs
		b <- rast(ext=as.vector(t(from@bbox)), nrow=from@grid@cells.dim[2], ncol=from@grid@cells.dim[1], crs=prj)
		if (inherits(from, "SpatialGridDataFrame")) {
			cls <- sapply(from@data, function(i) class(i)[1])
			if (all(cls == "numeric")) {
				nlyr(b) <- ncol(from@data)
				b <- setValues(b, as.matrix(from@data))
			} else {
				x <- vector(mode="list", length=ncol(from@data))
				for (i in 1:ncol(from@data)) {
					x[[i]] <- setValues(b, from@data[,i])
				}
				b <- rast(x)
			}
			names(b) <- colnames(from@data)
		} else {
			dim(b) <- c(from@grid@cells.dim[2], from@grid@cells.dim[1])
		}
		b
	}
)


setAs("SpatialPixels", "SpatRaster",
	function(from){
		if (methods::.hasSlot(from, "data")) {
			as(as(from, "SpatialGridDataFrame"), "SpatRaster")
		} else {
			as(as(from, "SpatialGrid"), "SpatRaster")
		}
	}
)



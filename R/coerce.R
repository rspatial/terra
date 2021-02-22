# Author: Robert J. Hijmans 
# Date : October 2018
# Version 1.0
# License GPL v3


#setMethod("as.list", signature(x="SpatRaster"), 
#	function(x, ...) {
#		lapply(1:nlyr(x), function(i) x[[i]])
#	}
#)
 

.as.image <- function(x, maxcells=10000) {
	x <- spatSample(x, size=maxcells, as.raster=TRUE)
	X <- xFromCol(x, 1:ncol(x))
	Y <- yFromRow(x, nrow(x):1)
	Z <- t(as.matrix(x, wide=TRUE)[nrow(x):1,]) 
	list(x=X, y=Y, z=Z)
}
 
 
setMethod("as.polygons", signature(x="SpatRaster"), 
	function(x, trunc=TRUE, dissolve=TRUE, values=TRUE, extent=FALSE) {
		p <- methods::new("SpatVector")
		if (extent) {
			p@ptr <- x@ptr$dense_extent()
		} else {
			opt <- spatOptions()
			p@ptr <- x@ptr$as_polygons(trunc[1], dissolve[1], values[1], TRUE, opt)
			#x <- messages(x)
		}
		messages(p, "as.polygons")
	}
)

setMethod("as.polygons", signature(x="SpatExtent"), 
	function(x, crs="") {
		p <- methods::new("SpatVector")
		p@ptr <- SpatVector$new(x@ptr, crs)
		messages(p, "as.polygons")
	}
)

setMethod("as.lines", signature(x="SpatExtent"), 
	function(x, crs="") {
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


setMethod("as.points", signature(x="SpatVector"), 
	function(x, multi=FALSE) {
		opt <- .getOptions()
		x@ptr <- x@ptr$as_points(multi)
		messages(x, "as.points")
	}
)


setMethod("as.points", signature(x="SpatRaster"), 
	function(x, values=TRUE) {
		p <- methods::new("SpatVector")
		opt <- .getOptions()
		p@ptr <- x@ptr$as_points(values, TRUE, opt)
		x <- messages(x, "as.points")
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


setMethod("as.matrix", signature(x="SpatRaster"), 
	function(x, wide=FALSE) {
		if (!hasValues(x)) {
			error("as.matrix", "SpatRaster has no cell values")
		}
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
)


setMethod("as.data.frame", signature(x="SpatRaster"), 
	function(x, xy=FALSE, cells=FALSE, na.rm=TRUE) {
		d <- NULL
		if (xy) {
			d <- xyFromCell(x, 1:ncell(x))
		} 
		if (cells) {
			d <- cbind(cell=1:ncell(x), d)
		}
		d <- cbind(d, values(x, mat=TRUE))
		if (na.rm) d <- stats::na.omit(d) 
		data.frame(d)
	}
)

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


# todo:
# for ncdf files (not yet natively supported in terra)
# check the variable to be used
# 
# check z values, other attributes such as NAvalue that may have been
# changed after creation of object from file
# RAT tables
# Author: Robert J. Hijmans 
# Date : February 2019
# Version 1.0
# License GPL v3

# todo:
# for ncdf files (not yet natively supported in terra)
# check the variable to be used
# 
# check z values, other attributes such as NAvalue that may have been
# changed after creation of object from file
# RAT tables

.fromRasterLayerBrick <- function(from) {
	f <- filename(from)
	if (f != "") {
		r <- rast(f)
		if (from@file@NAchanged) {
			warn("as,Raster", "changed NA value ignored")
		}
		return(r)
	} else {
		crsobj <- crs(from)
		if (is.na(crsobj)) {
			prj <- ""
		} else {
			crscom <- comment(crsobj)
			if (is.null(crscom)) {
				prj <- crsobj@projargs
			} else {
				prj <- crscom
			}
		}
		r <- rast(	nrows=nrow(from), 
					ncols=ncol(from),
					nlyrs=nlayers(from),
					crs=prj,
					extent=extent(from))
		if (hasValues(from)) {
			values(r) <- values(from)
		}
		names(r)  <- names(from)
	}
	return(r)
}

.fromRasterStack <- function(from) {
	x <- from[[1]]
	n <- nbands(x)
	if ((n > 1) & (n == nlayers(from))) {
		ff <- lapply(1:nlayers(from), function(i) { filename(from[[i]]) })
		if (length(unique(ff)) == 1) {
			r <- rast(filename(x))
			return(r)
		}
	} 
	s <- lapply(1:nlayers(from), function(i) {
		x <- from[[i]]
		.fromRasterLayerBrick(x)[[bandnr(x)]]
	})
	do.call(c, s)
}


setAs("Raster", "SpatRaster", 
	function(from) {
		if (inherits(from, "RasterLayer") | inherits(from, "RasterBrick")) { 
			.fromRasterLayerBrick(from)
		} else {
			.fromRasterStack(from)
		}
	}
)




setAs("SpatRaster", "Raster", 
	function(from) {
		s <- sources(from)
		nl <- nlyr(from)
		e <- as.vector(ext(from))
		prj <- .proj4(from)
		if (nl == 1) {
			if (s$source == "") {
				r <- raster(ncol=ncol(from), nrow=nrow(from), crs=prj,
			          xmn=e[1], xmx=e[2], ymn=e[3], ymx=e[4])
				if (hasValues(from)) {
					values(r) <- values(from)
				}
			} else {
				r <- raster(s$source)
			}
			names(r) <- names(from)
		} else {
			if (nrow(s) == 1 & s$source[1] != "") {
				r <- brick(s$source)
			} else if (all(s$source=="")) {
				r <- brick(ncol=ncol(from), nrow=nrow(from), crs=prj,
			          xmn=e[1], xmx=e[2], ymn=e[3], ymx=e[4], nl=nlyr(from))
				if (hasValues(from)) {
					values(r) <- values(from)
				}
			} else {
				x <- raster(ncol=ncol(from), nrow=nrow(from), crs=prj,
			          xmn=e[1], xmx=e[2], ymn=e[3], ymx=e[4])
				r <- list()
				for (i in 1:nl) {
					if (s$source[i] == "") {
						r[[i]] <- setValues(x, values(from[[i]]))
					} else {
						r[[i]] <- raster(s$source[i])
					}
				}
				r <- stack(r)
			}
		}
		return(r)
	}
)


# to sf from SpatVector
# could be in sf
.v2sf <- function(from) {
	txt <- 'sf::st_as_sf(as.data.frame(from, geom=TRUE), wkt="geometry", crs=from@ptr$get_crs("wkt"))'
	eval(parse(text = txt))
}

# from sf. first incomplete draft
.from_sf <- function(from) {
	sfi <- attr(from, "sf_column")
	geom <- from[[sfi]]
	crs <- attr(geom, "crs")$wkt
	attr(geom, "class") <- NULL
	types <- t(sapply(geom, function(i) attr(i, "class")))
	v <- list()
	for (i in 1:length(geom)) {
		if (inherits(geom[[i]], "POINT")) {
			v[[i]] <- cbind(i, 1, geom[[i]][1], geom[[i]][2], hole=0)				
		} else {
			vv <- list()
			for (j in 1:length(geom[[i]])) {
				if (inherits(geom[[i]][[j]], "list")) {
					vvv <- list()
					for (k in 1:length(geom[[i]][[j]])) {
						vvv[[k]] <- cbind(i, j, geom[[i]][[j]][[k]], hole=k-1) 
					}
					vv[[j]] <- do.call(rbind, vvv)
				} else {
					vv[[j]] <- cbind(i, j, geom[[i]][[j]], hole=0)
				}
			}
			v[[i]] <- do.call(rbind, vv)
		}
	}
	v <- do.call(rbind, v)
	if (ncol(v) > 5) {
		v <- cbind(v[,1:4], v[,ncol(v),drop=FALSE])
		warn("as", "Z/M dimension dropped")
	}
	colnames(v)[1:4] <- c("id", "part", "x", "y")
	types <- unique(gsub("MULTI", "", unique(types[,2])))
	if (length(types) > 1) {
		error("as,sf", "SpatVector currently only accepts one geometry type")
	}
	if (grepl("POINT", types, fixed=TRUE)) {
		gt = "points"
	} else if (grepl("LINE", types, fixed=TRUE)) {
		gt = "lines"
	} else if (grepl("POLY", types, fixed=TRUE)) {
		gt = "polygons"
	}
	if (ncol(from) > 1) {
		from[[sfi]] <- NULL
		d <- as.data.frame(from)
		vect(v, type=gt, att=d, crs=crs)
	} else {
		vect(v, type=gt, crs=crs)
	}
}


.from_sfc <- function(from) {
	geom = from
	v <- list()
	for (i in 1:length(geom)) {
		vv <- list()
		for (j in 1:length(geom[[i]])) {
			if (inherits(geom[[i]][[j]], "list")) {
				vvv <- list()
				for (k in 1:length(geom[[i]][[j]])) {
					vvv[[k]] <- cbind(i, j, geom[[i]][[j]][[k]], hole= k-1) 
				}
				vv[[j]] <- do.call(rbind, vvv)
			} else {
				vv[[j]] <- cbind(i, j, geom[[i]][[j]], hole=0) 
			}
		}
		v[[i]] <- do.call(rbind, vv)
	}
	v <- do.call(rbind, v)
	colnames(v)[1:4] <- c("id", "part", "x", "y")
	if (ncol(v) > 5) {
		v <- cbind(v[,1:4], v[,ncol(v),drop=FALSE])
		warn("as", "Z/M dimension dropped")
	}
	types <- class(from)[1]
	if (grepl("POINT", types, fixed=TRUE)) {
		gt = "points"
	} else if (grepl("LINE", types, fixed=TRUE)) {
		gt = "lines"
	} else if (grepl("POLY", types, fixed=TRUE)) {
		gt = "polygons"
	}
	vect(v, type=gt, crs="")
}



.from_sfg <- function(from) {
	geom = from
	v <- list()
	for (i in 1:length(geom)) {
		vv <- list()
		for (j in 1:length(geom[[i]])) {
			vv[[j]] <- cbind(i, j, geom[[i]][[j]], hole= j-1) 
		}
		v[[i]] <- do.call(rbind, vv)
	}
	v <- do.call(rbind, v)
	colnames(v)[1:4] <- c("id", "part", "x", "y")
	if (ncol(v) == 6) {
		v <- v[,-5]
	}
	types <- class(geom)[2]
	if (grepl("POINT", types, fixed=TRUE)) {
		gt = "points"
	} else if (grepl("LINE", types, fixed=TRUE)) {
		gt = "lines"
	} else if (grepl("POLY", types, fixed=TRUE)) {
		gt = "polygons"
	}
	vect(v, type=gt, crs="")
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

setAs("sfc", "SpatVector", 
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
		raster::geom(g, values(from), geomtype(from), crs(from))
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


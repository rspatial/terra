# Author: Robert J. Hijmans 
# Date : October 2018
# Version 1.0
# License GPL v3

 
setMethod("as.polygons", signature(x="SpatRaster"), 
	function(x, values=FALSE, na.rm=FALSE, ...) {
		p <- methods::new("SpatVector")
		p@ptr <- x@ptr$as_polygons(values, na.rm)
		x <- show_messages(x)
		show_messages(p)
	}
)

setMethod("as.lines", signature(x="SpatVector"), 
	function(x, ...) {
		x@ptr <- x@ptr$as_lines()
		show_messages(x)
	}
)

setMethod("as.points", signature(x="SpatRaster"), 
	function(x, values=FALSE, na.rm=FALSE, ...) {
		p <- methods::new("SpatVector")
		p@ptr <- x@ptr$as_points(values, na.rm)
		x <- show_messages(x)
		show_messages(p)
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
	function(x, ...) {
		paste( x@ptr$vector, collapse=", ")
	}
)


setMethod("as.vector", signature(x="SpatRaster"), 
	function(x, mode="any") {
		values(x, FALSE)
	}
)


setMethod("as.matrix", signature(x="SpatRaster"), 
	function(x, wide=FALSE, ...) {
		if (!hasValues(x)) {
			stop("SpatRaster has no cell values")
		}
		if (wide) {
			if (nlyr(x) > 1) {
				m <- values(x, matrix=TRUE)
				m <- lapply(1:ncol(m), function(i) {
					matrix(m[,i],nrow=nrow(x),byrow=TRUE)
					})
				m <- do.call(cbind, m)
			} else {
				m <- matrix(values(x, matrix=FALSE),nrow=nrow(x),byrow=TRUE)
			}
		} else {
			m <- values(x, matrix=TRUE)
		}
		m
	}
)


setMethod("as.data.frame", signature(x="SpatRaster"), 
	function(x, xy=FALSE, cells=FALSE, ...) {
		d <- NULL
		if (xy) {
			d <- xyFromCell(x, 1:ncell(x))
		} 
		if (cells) {
			d <- cbind(cell=1:ncell(x), d)
		}
		d <- cbind(d, values(x, matrix=TRUE))
		data.frame(d)
	}
)



setMethod("as.array", signature(x="SpatRaster"), 
	function(x, ...) {
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

.RasterLayerToSpatRaster <- function(from) { 

	f <- filename(from)
	if (f != "") {
		# a bit more tricky with ncdf...
		r <- rast(f)
	} else {
		e <- extent(from)
		r <- rast(ncol=ncol(from), nrow=nrow(from), crs=crs(from),
		          xmin=e@xmin, xmax=e@xmax, ymin=e@ymin, ymax=e@ymax,
				  nlyr=1)				
		if (hasValues(from)) {
			values(r) <- values(from)
		}
	}
	names(r) <- names(from)
	return(r)
}



.RasterBrickToSpatRaster <- function(from) { 

	f <- filename(from)
	if (f != "") {
		# a bit more tricky with ncdf...
		if (from@file@NAchanged) {
			warning("changed NAvalue ignored")
		}
		r <- rast(f)
	} else {
		e <- extent(from)
		nl <- nlayers(from)
		r <- rast(ncol=ncol(from), nrow=nrow(from), crs=crs(from),
		          xmin=e@xmin, xmax=e@xmax, ymin=e@ymin, ymax=e@ymax,
				  nlyr=nl)
		if (hasValues(from)) {
			values(r) <- values(from)
		}
	}
	names(r) <- names(from)
	return(r)
}



.RasterStackToSpatRaster <- function(from) { 
	nl <- nlayers(from)
	rr <- methods::as("SpatRaster", from[[1]])		
	nb <- nbands(from)
	
	if ((nb > 1) & (nb == nl)) {
		ff <- lapply(1:nlayers(from), function(i) { filename(from[[i]]) })
		ff <- unique(ff)
		if ((length(ff) == 1) & (ff != "")) {
			rr <- rast(ff)
		}
		names(rr) <- names(from)
		return(rr)
	}
	
	if (nl > 1) {
		for (i in 2:nl) {
			rr <- c(rr, methods::as("SpatRaster", from[[i]]))
		}
	}
	names(rr) <- names(from)	
	rr
}


setAs("Raster", "SpatRaster", 
	function(from) { 
		if (inherits(from, "RasterLayer")) {
			r <- .RasterLayerToSpatRaster(from)
		}
		if (inherits(from, "RasterStack")) {
			r <- .RasterStackToSpatRaster(from)
		}
		if (inherits(from, "RasterBrick")) {
			r <- .RasterBrickToSpatRaster(from)
		}
		show_messages(r, "coerce")
	}
)

setAs("SpatRaster", "Raster", 
	function(from) {
		s <- sources(from)
		nl <- nlyr(from)
		e <- as.vector(ext(from))
		if (nl == 1) {
			if (s$source == "") {
				r <- raster(ncol=ncol(from), nrow=nrow(from), crs=crs(from),
			          xmn=e[1], xmx=e[2], ymn=e[3], ymx=e[4])
				if (.hasValues(from)) {
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
				r <- brick(ncol=ncol(from), nrow=nrow(from), crs=crs(from),
			          xmn=e[1], xmx=e[2], ymn=e[3], ymx=e[4], nl=nlyr(from))
				if (.hasValues(from)) {
					values(r) <- values(from)
				}
			} else {
				x <- raster(ncol=ncol(from), nrow=nrow(from), crs=crs(from),
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



setAs("SpatVector", "Spatial", 
	function(from) {
		g <- geom(from)
		colnames(g)[1] <- "object"
		raster::geom(g, values(from), geomtype(from), crs(from))
	}
)


setAs("Spatial", "SpatVector", 
	function(from) {
		g <- geom(from)
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
			g <- cbind(g[,1,drop=FALSE], part=1:nrow(g), g[,2:3])
		}
		if (methods::.hasSlot(from, "data")) {
			v <- vect(g, vtype, from@data, crs(from))
		} else {
			v <- vect(g, vtype, crs=crs(from))
		}
		return(v)
	}
)

